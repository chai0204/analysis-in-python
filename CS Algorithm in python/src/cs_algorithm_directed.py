# cs_algorithm_directed.py
"""
目的：
CSアルゴリズム(Isozaki, 2014)で提案されているロジックに準拠し、データから部分的有向非巡回グラフ（PDAG）を構築します。
骨格発見、V構造特定、矛盾解決、信頼できない向きの削除、論理ルールによる伝播の全ステップを実行します。
最終的に、各辺の強さをバックドア基準（親ノード）を統制変数とした偏相関係数として算出します。

設定項目：
def main() 内の設定項目を編集して実行してください。

入力：
- CSVファイル

出力：
- コンソールへの分析過程と最終結果の表示
- 分析結果のJSONファイル
"""

import pandas as pd
import pingouin as pg
from itertools import combinations
from collections import defaultdict
import networkx as nx
import traceback
import json

# --- ヘルパー関数 ---

def get_v_structure_tuples(directed_edges: set):
    """V構造の有向辺セットを、可読な(X, Y, Z)のタプルセットに変換する"""
    colliders = defaultdict(list)
    for u, v in directed_edges: colliders[v].append(u)
    v_structures = set()
    for z, parents in colliders.items():
        if len(parents) >= 2:
            for x, y in combinations(parents, 2):
                # BUG FIX: list [z] ではなく tuple (z,) を結合する
                v_structures.add(tuple(sorted((x,y))) + (z,))
    return v_structures

# --- フェーズ1：骨格発見 ---

def discover_skeleton(df: pd.DataFrame, alpha: float, max_control_vars: int):
    """CSアルゴリズムに基づき、グラフの骨格と、向き付けに必要な分離集合・p値を発見する"""
    print("\n--- [フェーズ1] グラフ骨格の発見 ---")
    variables = list(df.columns)
    G = nx.complete_graph(variables)
    sepsets = defaultdict(list)
    sepset_pvals = {}

    initial_edges = G.number_of_edges()
    print(f"  - 分析開始時のグラフ: 完全グラフ (辺の数: {initial_edges})")

    print("\n[ステップ1.1] 0次の独立性検定")
    edges_before = G.number_of_edges()
    pairwise_results = pg.pairwise_corr(df, method='pearson')
    for _, row in pairwise_results.iterrows():
        x, y, p_val = row['X'], row['Y'], row['p-unc']
        if p_val > alpha:
            if G.has_edge(x,y): 
                G.remove_edge(x, y)
                print(f"  - [辺の削除] {x} - {y} (p={p_val:.4f})")
            key = tuple(sorted((x, y)))
            sepsets[key].append([])
            sepset_pvals[key] = p_val
    edges_after = G.number_of_edges()
    print(f"  - [結果] 削除された辺の数: {edges_before - edges_after} | 残りの辺の数: {edges_after}")

    temp_v_structures = find_v_structures(G, sepsets)
    for n in range(1, max_control_vars + 1):
        print(f"\n[ステップ1.2] {n}次の条件付き独立性検定")
        edges_before_n = G.number_of_edges()

        if temp_v_structures:
            print(f"  - 現在のV構造（MBCチェック用）: { {f'{u}->{z}<-{v}' for u,z,v in get_v_structure_tuples(temp_v_structures)} }")

        edges_removed_in_this_round = False
        for x, y in sorted(list(G.edges())):
            potential_S = (set(G.neighbors(x)) | set(G.neighbors(y))) - {x, y}
            if len(potential_S) < n: continue
            for s in combinations(potential_S, n):
                if check_strict_mbc(G, x, y, s, temp_v_structures):
                    continue
                
                pcorr = pg.partial_corr(data=df, x=x, y=y, covar=list(s))
                p_val = pcorr['p-val'].iloc[0]
                if p_val > alpha:
                    print(f"  - [辺の削除] {x} - {y} | {s} (p={p_val:.4f})")
                    if G.has_edge(x,y): G.remove_edge(x, y)
                    key = tuple(sorted((x, y))); sepsets[key].append(list(s)); sepset_pvals[key] = p_val
                    edges_removed_in_this_round = True
                    break
        
        edges_after_n = G.number_of_edges()
        print(f"  - [結果] このステップで削除された辺の数: {edges_before_n - edges_after_n} | 残りの辺の数: {edges_after_n}")

        if edges_removed_in_this_round:
            temp_v_structures = find_v_structures(G, sepsets)
        else:
            print("  - 辺の削除がなかったため、骨格発見を完了します。")
            break
    print(f"\n--- 骨格発見 完了（最終的な辺の数: {G.number_of_edges()}） ---")
    return G, sepsets, sepset_pvals

def find_v_structures(G: nx.Graph, sepsets: dict):
    directed_edges = set()
    variables = list(G.nodes())
    for x, y in combinations(variables, 2):
        if not G.has_edge(x, y):
            common_neighbors = set(G.neighbors(x)) & set(G.neighbors(y))
            for z in common_neighbors:
                if not any(z in s for s in sepsets.get(tuple(sorted((x,y))), [])):
                    directed_edges.add((x, z)); directed_edges.add((y, z))
    return directed_edges

def check_strict_mbc(G: nx.Graph, x: str, y: str, s: tuple, directed_edges: set):
    for z in s:
        if not G.has_edge(x, z) or not G.has_edge(y, z): continue
        if not ((x, z) in directed_edges and (y, z) in directed_edges): return False
    return True

# --- フェーズ2：向き付け ---

def orient_graph(G: nx.Graph, sepsets: dict, sepset_pvals: dict, df: pd.DataFrame, alpha: float):
    """骨格グラフに対し、向き付けのルールを適用してPDAG（部分的有向非巡回グラフ）を返す"""
    print("\n--- [フェーズ2] エッジの向き付け ---")
    
    print("\n[ステップ2.1] V構造の特定")
    directed_edges = find_v_structures(G, sepsets)
    v_tuples = get_v_structure_tuples(directed_edges)
    if v_tuples: print(f"  - 発見されたV構造: { {f'{u}->{z}<-{v}' for u,z,v in v_tuples} }")
    else: print("  - V構造は見つかりませんでした。")

    directed_edges = resolve_inconsistencies(directed_edges, G, sepset_pvals)
    directed_edges = handle_unreliable_directions(directed_edges, G, df, alpha)
    directed_edges = apply_orientation_rules(G, directed_edges)
    
    final_undirected = {tuple(sorted(e)) for e in G.edges()}
    final_directed = set()
    processed_bidirectional = set()
    for u, v in directed_edges.copy():
        key = tuple(sorted((u,v)))
        if (v, u) in directed_edges:
            if key not in processed_bidirectional: processed_bidirectional.add(key)
        else:
            final_directed.add((u,v))
            if key in final_undirected: final_undirected.remove(key)
    final_undirected.update(processed_bidirectional)
    final_directed = {(u,v) for u,v in final_directed if tuple(sorted((u,v))) not in processed_bidirectional}
    print("\n--- 向き付け 完了 ---")
    return final_directed, final_undirected

def resolve_inconsistencies(directed_edges: set, G: nx.Graph, sepset_pvals: dict):
    """
    矛盾する双方向エッジを、サイクルを生成しないように解決する。
    サイクルが生成される場合は、より安全な無向化を選択する。
    """
    bi_directional_pairs = {tuple(sorted((u, v))) for u, v in directed_edges if (v, u) in directed_edges}
    if not bi_directional_pairs:
        return directed_edges
    
    print("\n[ステップ2.2] 矛盾の解決（サイクルチェック実行）")
    for u, v in sorted(list(bi_directional_pairs)):
        # このペアがまだ処理対象か再チェック（ループ内で集合が変更されるため）
        if not ((u, v) in directed_edges and (v, u) in directed_edges):
            continue

        print(f"  - [矛盾検出] {u} <--> {v}")
        
        # --- p値の収集ロジックは変更なし ---
        p_val_uv, cause_uv = 0, None
        for w in G.neighbors(v):
            if w != u and (w, v) in directed_edges and not G.has_edge(w, u):
                key = tuple(sorted((w, u)))
                if key in sepset_pvals and sepset_pvals[key] > p_val_uv:
                    p_val_uv = sepset_pvals[key]; cause_uv = f"{w} -> {v} <- {u}"
        
        p_val_vu, cause_vu = 0, None
        for z in G.neighbors(u):
            if z != v and (z, u) in directed_edges and not G.has_edge(z, v):
                key = tuple(sorted((z, v)))
                if key in sepset_pvals and sepset_pvals[key] > p_val_vu:
                    p_val_vu = sepset_pvals[key]; cause_vu = f"{z} -> {u} <- {v}"

        print(f"    - u -> v の根拠: {cause_uv or '不明'} (p値: {p_val_uv:.4f})")
        print(f"    - v -> u の根拠: {cause_vu or '不明'} (p値: {p_val_vu:.4f})")

        # --- サイクルチェックを組み込んだ解決ロジック ---
        base_edges = directed_edges - {(u, v), (v, u)}
        temp_graph = nx.DiGraph()
        temp_graph.add_nodes_from(G.nodes())
        temp_graph.add_edges_from(list(base_edges))

        # u -> v を残す場合 (p_val_uv > p_val_vu)
        if p_val_uv > p_val_vu:
            # v から u へのパスが既に存在するか？ (存在すれば u->v を追加するとサイクルになる)
            if nx.has_path(temp_graph, v, u):
                print(f"    - [解決] {u} -> {v} はサイクルを生成するため、両方の向きを削除（無向化）")
                directed_edges.remove((u, v))
                directed_edges.remove((v, u))
            else:
                print(f"    - [解決] {v} -> {u} を削除")
                directed_edges.remove((v, u))
        
        # v -> u を残す場合 (p_val_vu > p_val_uv)
        elif p_val_vu > p_val_uv:
            # u から v へのパスが既に存在するか？ (存在すれば v->u を追加するとサイクルになる)
            if nx.has_path(temp_graph, u, v):
                print(f"    - [解決] {v} -> {u} はサイクルを生成するため、両方の向きを削除（無向化）")
                directed_edges.remove((u, v))
                directed_edges.remove((v, u))
            else:
                print(f"    - [解決] {u} -> {v} を削除")
                directed_edges.remove((u, v))
        
        # p値が同等または不明な場合
        else:
            print("    - [解決] p値が同等または不明なため、両方の向きを削除（無向化）")
            directed_edges.remove((u, v))
            directed_edges.remove((v, u))
            
    return directed_edges

def handle_unreliable_directions(directed_edges: set, G: nx.Graph, df: pd.DataFrame, alpha: float):
    colliders = defaultdict(list); [colliders[v].append(u) for u, v in directed_edges]
    found = False
    for z, parents in sorted(colliders.items()):
        if len(parents) < 2: continue
        for x, y in combinations(sorted(parents), 2):
            for w, w_parents in sorted(colliders.items()):
                if w == z or x not in w_parents or y not in w_parents: continue
                pcorr = pg.partial_corr(data=df, x=x, y=y, covar=[z])
                if pcorr is not None and not pcorr.empty and pcorr['p-val'].iloc[0] > alpha:
                    if not found: print("\n[ステップ2.3] 信頼できない向きの処理")
                    print(f"  - [パターン発見] {x}->{z}<-{y} と {x}->{w}<-{y}")
                    print(f"    (理由: {x}と{y}が{z}で条件付き独立 p={pcorr['p-val'].iloc[0]:.4f})")
                    print(f"    - [修正] {x}->{z} と {y}->{z} の向きを削除")
                    if (x, z) in directed_edges: directed_edges.remove((x, z))
                    if (y, z) in directed_edges: directed_edges.remove((y, z))
                    found = True
    return directed_edges

def apply_orientation_rules(G: nx.Graph, directed_edges: set):
    """
    論理ルール(R1-R4)に基づき、サイクルを生成しないように向き付けを伝播させる。
    Meek, C. (1995) Causal inference from graphical models.
    """
    print("\n[ステップ2.4] 論理ルールに基づく向き付けの伝播（R1-R4, サイクルチェック実行）")
    
    def _has_path(source, target, edges):
        temp_graph = nx.DiGraph()
        temp_graph.add_nodes_from(G.nodes())
        temp_graph.add_edges_from(list(edges))
        return nx.has_path(temp_graph, source, target)

    while True:
        new_orientations_found = False
        
        # --- ルール1 (R1): X -> Y - Z (X,Zが非隣接) => Y -> Z ---
        for x, y in sorted(list(directed_edges)):
            if (y, x) in directed_edges: continue
            for z in sorted(list(G.neighbors(y))):
                if z != x and not G.has_edge(x, z) and (y, z) not in directed_edges and (z, y) not in directed_edges:
                    if not _has_path(z, y, directed_edges):
                        print(f"  - [ルール1適用] {x} -> {y} - {z} (かつ {x},{z}は非隣接) => {y} -> {z}")
                        directed_edges.add((y, z))
                        new_orientations_found = True
                    else:
                        print(f"  - [ルール1 スキップ] {y} -> {z} はサイクルを生成するため適用しません")

        # --- ルール2 (R2): X -> Y -> Z (X-Z) => X -> Z ---
        for x, y in sorted(list(directed_edges)):
            if (y, x) in directed_edges: continue
            for z in sorted(list(G.neighbors(y))):
                if z != x and (y, z) in directed_edges and (z, y) not in directed_edges and \
                   G.has_edge(x, z) and (x, z) not in directed_edges and (z, x) not in directed_edges:
                    if not _has_path(z, x, directed_edges):
                        print(f"  - [ルール2適用] {x} -> {y} -> {z} (かつ {x}-{z}) => {x} -> {z}")
                        directed_edges.add((x, z))
                        new_orientations_found = True
                    else:
                        print(f"  - [ルール2 スキップ] {x} -> {z} はサイクルを生成するため適用しません")

        # --- ルール3 (R3): X-Y, X-Z, Y->W, Z->W (Y,Z非隣接) => X->W ---
        # Y->W<-Z というV構造を探す
        colliders = defaultdict(list)
        for u, v in directed_edges:
            if (v, u) not in directed_edges:
                colliders[v].append(u)
        
        for w, parents in sorted(colliders.items()):
            if len(parents) < 2: continue
            for y, z in combinations(sorted(parents), 2):
                # 親同士(Y,Z)が非隣接かチェック
                if not G.has_edge(y, z):
                    # Y,Zに共通の隣接ノードXを探す
                    common_neighbors_of_yz = set(G.neighbors(y)) & set(G.neighbors(z))
                    for x in sorted(list(common_neighbors_of_yz)):
                        if x != w and (x, w) not in directed_edges and (w, x) not in directed_edges:
                            if not _has_path(w, x, directed_edges):
                                print(f"  - [ルール3適用] {y}->{w}<-{z} と {y}-{x}-{z} => {x} -> {w}")
                                directed_edges.add((x, w))
                                new_orientations_found = True
                            else:
                                print(f"  - [ルール3 スキップ] {x} -> {w} はサイクルを生成するため適用しません")

        # --- ルール4 (R4): X->Y->Z, X-W-Z => W->Z ---
        # X->Y->Z パスを探す
        for x, y in sorted(list(directed_edges)):
            if (y, x) in directed_edges: continue
            for z_node in sorted(list(G.neighbors(y))):
                if (y, z_node) in directed_edges and (z_node, y) not in directed_edges and x != z_node:
                    # X-W-Z パスを探す
                    common_neighbors_of_xz = (set(G.neighbors(x)) & set(G.neighbors(z_node))) - {y}
                    for w in sorted(list(common_neighbors_of_xz)):
                        if (w, z_node) not in directed_edges and (z_node, w) not in directed_edges:
                            if not _has_path(z_node, w, directed_edges):
                                print(f"  - [ルール4適用] {x}->{y}->{z_node} と {x}-{w}-{z_node} => {w} -> {z_node}")
                                directed_edges.add((w, z_node))
                                new_orientations_found = True
                            else:
                                print(f"  - [ルール4 スキップ] {w} -> {z_node} はサイクルを生成するため適用しません")

        if not new_orientations_found:
            break
            
    return directed_edges


# --- フェーズ3：強さ計算と結果表示 ---

def calculate_and_summarize(df: pd.DataFrame, directed_edges: set, undirected_edges: set, alpha: float, output_json_path: str):
    """有向グラフの各辺に対し、バックドア基準で偏相関係数を計算し、結果を要約・JSON出力する"""
    print("\n--- [フェーズ3] パスの強さの計算と最終サマリー ---")
    print("\n[ステップ3.1] パスの強さの計算（バックドア基準）")
    parents = defaultdict(set); [parents[v].add(u) for u, v in directed_edges]
    final_strengths = []
    edges_to_process = sorted(list(directed_edges)) + sorted(list(undirected_edges))
    for u, v in edges_to_process:
        control_vars = list((parents[u] | parents[v]) - {u, v})
        try:
            if not control_vars:
                res = pg.corr(df[u], df[v]); strength = res['r'].iloc[0]; p_val = res['p-val'].iloc[0]
                print(f"  - {u} -- {v}: 相関係数 = {strength:.3f} (p={p_val:.4f})")
            else:
                res = pg.partial_corr(data=df, x=u, y=v, covar=control_vars); strength = res['r'].iloc[0]; p_val = res['p-val'].iloc[0]
                print(f"  - {u} -- {v}: 偏相関係数 = {strength:.3f} (p={p_val:.4f}), 統制変数: {control_vars}")
            final_strengths.append({'edge': (u, v), 'strength': strength, 'p_value': p_val, 'controls': control_vars})
        except Exception as e:
            print(f"  - {u} --- {v} の計算でエラー: {e}")

    print("\n\n--- ★★★ 分析結果の最終サマリー ★★★ ---")
    results_map = {res['edge']: res for res in final_strengths if res['p_value'] < alpha}
    if not results_map: 
        print("\n統計的に有意なパスは見つかりませんでした。")
        return
    
    bi, uni, undir = set(), set(), set(); processed_bi = set()
    json_output = []

    for u, v in directed_edges:
        if (u, v) in results_map:
            if (v, u) in directed_edges:
                key = tuple(sorted((u, v)))
                if key not in processed_bi: bi.add(key); processed_bi.add(key)
            else: uni.add((u, v))
    for u, v in undirected_edges:
        if (u, v) in results_map:
            key = tuple(sorted((u, v)))
            if key not in processed_bi: undir.add(key)

    print(f"\n【合計: {len(bi) + len(uni) + len(undir)}本の有意なパスが特定されました】")
    # 双方向パス
    print(f"\n[双方向パス: {len(bi)}本]")
    if bi:
        for u, v in sorted(list(bi)):
            res = results_map.get((u, v)) or results_map.get((v, u))
            if res: 
                print(f"  - {u} <--> {v} (強さ: {res['strength']:.3f}, p値: {res['p_value']:.4f}, 統制変数: {res['controls'] if res['controls'] else 'なし'})")
                json_output.append({
                    "変数1": u, "向き": "<-->", "変数2": v,
                    "偏相関係数": res['strength'], "p値": res['p_value'], "統制変数群": res['controls']
                })
    else: print("  - なし")
    # 単方向パス
    print(f"\n[単方向パス: {len(uni)}本]")
    if uni:
        for u, v in sorted(list(uni)):
            res = results_map.get((u, v))
            if res: 
                print(f"  - {u} --> {v} (強さ: {res['strength']:.3f}, p値: {res['p_value']:.4f}, 統制変数: {res['controls'] if res['controls'] else 'なし'})")
                json_output.append({
                    "変数1": u, "向き": "-->", "変数2": v,
                    "偏相関係数": res['strength'], "p値": res['p_value'], "統制変数群": res['controls']
                })
    else: print("  - なし")
    # 方向未決定パス
    print(f"\n[有意だが方向未決定のパス: {len(undir)}本]")
    if undir:
        for u, v in sorted(list(undir)):
            res = results_map.get((u, v))
            if res: 
                print(f"  - {u} --- {v} (強さ: {res['strength']:.3f}, p値: {res['p_value']:.4f}, 統制変数: {res['controls'] if res['controls'] else 'なし'})")
                json_output.append({
                    "変数1": u, "向き": "---", "変数2": v,
                    "偏相関係数": res['strength'], "p値": res['p_value'], "統制変数群": res['controls']
                })
    else: print("  - なし")

    # --- JSONファイルへの書き出し ---
    if output_json_path:
        try:
            with open(output_json_path, 'w', encoding='utf-8') as f:
                json.dump(json_output, f, ensure_ascii=False, indent=2)
            print(f"\n分析結果が '{output_json_path}' にJSON形式で保存されました。")
        except Exception as e:
            print(f"\nJSONファイルへの書き出し中にエラーが発生しました: {e}")


# --- 実行ブロック ---

def run_directed_analysis(input_csv_path: str, significance_level: float = 0.05, max_control_vars: int = 4, output_json_path: str = None):
    """
    有向グラフ分析を実行するメイン関数。

    Args:
        input_csv_path (str): 分析対象データ（CSV形式）のファイルパス。
        significance_level (float, optional): 統計的検定の有意水準（α）。デフォルトは 0.05。
        max_control_vars (int, optional): 条件付き独立性検定で考慮する最大変数数。デフォルトは 4。
        output_json_path (str, optional): 結果をJSON形式で保存するファイルパス。デフォルトは None。
    """
    try:
        df = pd.read_csv(input_csv_path, encoding='utf-8')
        print(f"CSVファイル '{input_csv_path}' の読み込みに成功しました。")
        
        # フェーズ1: 骨格発見
        G, sepsets, sepset_pvals = discover_skeleton(df, significance_level, max_control_vars)

        # フェーズ2: 向き付け
        final_directed, final_undirected = orient_graph(G, sepsets, sepset_pvals, df, significance_level)

        # フェーズ3: 強さ計算と結果表示
        calculate_and_summarize(df, final_directed, final_undirected, significance_level, output_json_path)

    except FileNotFoundError:
        print(f"エラー: ファイル '{input_csv_path}' が見つかりません。パスを確認してください。")
    except Exception as e:
        print("予期せぬエラーが発生しました。")
        traceback.print_exc()

def main():
    """スクリプトのメイン処理（設定項目を編集して実行）"""
    # --- 設定項目 ---
    INPUT_CSV_PATH = 'data/sample_data.csv' #任意の入力データのパスを入力
    SIGNIFICANCE_LEVEL = 0.05 #有意水準α
    MAX_CONTROL_VARS = 4 #条件付き独立性検定で考慮する最大変数数
    OUTPUT_JSON_PATH = 'output/causal_analysis_results.json' # 出力ファイル名/パス

    # 分析実行
    run_directed_analysis(
        input_csv_path=INPUT_CSV_PATH,
        significance_level=SIGNIFICANCE_LEVEL,
        max_control_vars=MAX_CONTROL_VARS,
        output_json_path=OUTPUT_JSON_PATH
    )

if __name__ == '__main__':
    main()
