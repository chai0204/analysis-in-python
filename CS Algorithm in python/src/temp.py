# temp.py
"""
目的：
CSアルゴリズム(Isozaki, 2014)に基づき、データから無向グラフを構築します。
骨格発見のステップのみを実行し、マルコフブランケットに基づき各辺の強さを偏相関係数として算出します。

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
                v_structures.add(tuple(sorted((x,y))) + (z,))
    return v_structures

# --- フェーズ1：骨格発見 ---

def discover_skeleton(df: pd.DataFrame, alpha: float, max_control_vars: int):
    """CSアルゴリズムに基づき、グラフの骨格を発見する"""
    print("\n--- [フェーズ1] グラフ骨格の発見 ---")
    variables = list(df.columns)
    G = nx.complete_graph(variables)
    sepsets = defaultdict(list)

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
    edges_after = G.number_of_edges()
    print(f"  - [結果] 削除された辺の数: {edges_before - edges_after} | 残りの辺の数: {edges_after}")

    temp_v_structures = find_v_structures(G, sepsets)
    for n in range(1, max_control_vars + 1):
        print(f"\n[ステップ1.2] {n}次の条件付き独立性検定")
        edges_before_n = G.number_of_edges()

        if temp_v_structures:
            v_tuples = get_v_structure_tuples(temp_v_structures)
            if v_tuples:
                print(f"  - 現在のV構造（MBCチェック用）: { {f'{u}->{z}<-{v}' for u,v,z in v_tuples} }")


        edges_removed_in_this_round = False
        edges_to_check = sorted(list(G.edges()))
        for x, y in edges_to_check:
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
                    key = tuple(sorted((x, y))); sepsets[key].append(list(s))
                    edges_removed_in_this_round = True
                    break 
            if not G.has_edge(x, y):
                continue

        edges_after_n = G.number_of_edges()
        print(f"  - [結果] このステップで削除された辺の数: {edges_before_n - edges_after_n} | 残りの辺の数: {edges_after_n}")

        if edges_removed_in_this_round:
            temp_v_structures = find_v_structures(G, sepsets)
        else:
            print("  - 辺の削除がなかったため、骨格発見を完了します。")
            break
    print(f"\n--- 骨格発見 完了（最終的な辺の数: {G.number_of_edges()}） ---")
    return G

def find_v_structures(G: nx.Graph, sepsets: dict):
    """骨格発見のMBCチェックで使うための一時的なV構造発見"""
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
    """骨格発見のMBCチェック"""
    for z in s:
        if not G.has_edge(x, z) or not G.has_edge(y, z): continue
        if not ((x, z) in directed_edges and (y, z) in directed_edges): return False
    return True

# --- フェーズ2：強さ計算と結果表示 (無向グラフ用) ---

def calculate_and_summarize_undirected(df: pd.DataFrame, G: nx.Graph, alpha: float, output_json_path: str):
    """
    無向グラフの各辺に対し、隣接ノードを統制変数として偏相関係数を計算し、結果を要約・JSON出力する。
    これは有向グラフ版のバックドア基準のアナロジーです。
    """
    print("\n--- [フェーズ2] パスの強さの計算と最終サマリー ---")
    print("\n[ステップ2.1] パスの強さの計算（隣接ノード基準）")
    
    final_strengths = []
    undirected_edges = sorted([tuple(sorted(e)) for e in G.edges()])

    for u, v in undirected_edges:
        # 無向グラフにおけるマルコフブランケット(隣接ノード)の考え方を応用し、
        # uとv両方のマルコフブランケットの和集合を統制変数とする。
        # これにより、u-v間の交絡となりうるパスの影響を最大限除去する。
        control_vars = list((set(G.neighbors(u)) | set(G.neighbors(v))) - {u, v})
        
        try:
            if not control_vars:
                res = pg.corr(df[u], df[v])
                strength = res['r'].iloc[0]
                p_val = res['p-val'].iloc[0]
                print(f"  - {u} -- {v}: 相関係数 = {strength:.3f} (p={p_val:.4f})")
            else:
                res = pg.partial_corr(data=df, x=u, y=v, covar=control_vars)
                strength = res['r'].iloc[0]
                p_val = res['p-val'].iloc[0]
                print(f"  - {u} -- {v}: 偏相関係数 = {strength:.3f} (p={p_val:.4f}), 統制変数: {control_vars}")
            
            final_strengths.append({'edge': (u, v), 'strength': strength, 'p_value': p_val, 'controls': control_vars})
        except Exception as e:
            # 統制変数が多すぎる場合（multicollinearity等）のエラーハンドリング
            print(f"  - {u} --- {v} の計算でエラー: {e}。統制変数なしで再計算します。")
            try:
                res = pg.corr(df[u], df[v])
                strength = res['r'].iloc[0]
                p_val = res['p-val'].iloc[0]
                print(f"    - {u} -- {v}: 相関係数 = {strength:.3f} (p={p_val:.4f})")
                final_strengths.append({'edge': (u, v), 'strength': strength, 'p_value': p_val, 'controls': []})
            except Exception as e2:
                print(f"    - 再計算でもエラー: {e2}")

    print("\n\n--- ★★★ 分析結果の最終サマリー ★★★ ---")
    
    significant_edges = [res for res in final_strengths if res['p_value'] < alpha]

    if not significant_edges:
        print("\n統計的に有意なパスは見つかりませんでした。")
        if output_json_path:
            try:
                with open(output_json_path, 'w', encoding='utf-8') as f:
                    json.dump([], f, ensure_ascii=False, indent=2)
                print(f"\n空の結果が '{output_json_path}' にJSON形式で保存されました。")
            except Exception as e:
                print(f"\nJSONファイルへの書き出し中にエラーが発生しました: {e}")
        return

    json_output = []
    print(f"\n【合計: {len(significant_edges)}本の有意なパスが特定されました】")
    for res in sorted(significant_edges, key=lambda x: x['edge']):
        u, v = res['edge']
        print(f"  - {u} --- {v} (強さ: {res['strength']:.3f}, p値: {res['p_value']:.4f}, 統制変数: {res['controls'] if res['controls'] else 'なし'})")
        json_output.append({
            "変数1": u, "向き": "---", "変数2": v,
            "偏相関係数": res['strength'], "p値": res['p_value'], "統制変数群": res['controls']
        })

    # --- JSONファイルへの書き出し ---
    if output_json_path:
        try:
            with open(output_json_path, 'w', encoding='utf-8') as f:
                json.dump(json_output, f, ensure_ascii=False, indent=2)
            print(f"\n分析結果が '{output_json_path}' にJSON形式で保存されました。")
        except Exception as e:
            print(f"\nJSONファイルへの書き出し中にエラーが発生しました: {e}")


# --- 実行ブロック ---

def run_undirected_analysis(input_csv_path: str, significance_level: float = 0.05, max_control_vars: int = 4, output_json_path: str = None):
    """
    無向グラフ分析を実行するメイン関数。

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
        G = discover_skeleton(df, significance_level, max_control_vars)

        # フェーズ2: 強さ計算と結果表示
        calculate_and_summarize_undirected(df, G, significance_level, output_json_path)

    except FileNotFoundError:
        print(f"エラー: ファイル '{input_csv_path}' が見つかりません。パスを確認してください。")
    except Exception as e:
        print("予期せぬエラーが発生しました。")
        traceback.print_exc()

def main():
    """スクリプトのメイン処理（設定項目を編集して実行）"""
    # --- 設定項目 ---
    INPUT_CSV_PATH = 'data/cleaned_no_outliers3.csv' # or 'data/sachs_data.csv'
    SIGNIFICANCE_LEVEL = 0.05
    MAX_CONTROL_VARS = 4
    OUTPUT_JSON_PATH = 'output/causal_analysis_results_undirected.json' # 出力ファイル名

    # 分析実行
    run_undirected_analysis(
        input_csv_path=INPUT_CSV_PATH,
        significance_level=SIGNIFICANCE_LEVEL,
        max_control_vars=MAX_CONTROL_VARS,
        output_json_path=OUTPUT_JSON_PATH
    )

if __name__ == '__main__':
    main()
