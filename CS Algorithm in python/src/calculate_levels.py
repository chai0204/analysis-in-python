
import networkx as nx
from collections import defaultdict

# 分析結果の有向パスデータ
directed_edges = [
    ("CFC", "MW_S"), ("CFC", "TIPI_Ne"), ("MRS_co", "MRS_av"),
    ("MRS_co", "MRS_ef"), ("MRS_co", "MRS_in"), ("MRS_ef", "Po"),
    ("MRS_in", "CFA"), ("MRS_in", "SSCS_CPI"), ("MRS_re", "MRS_co"),
    ("MRS_re", "MRS_ef"), ("MRS_re", "Te"), ("MW_S", "MW_D"),
    ("MW_S", "TIPI_Ex"), ("Po", "TIPI_Op"), ("SSCS_CPI", "SSCS_CSE"),
    ("SSCS_CSE", "CFC"), ("SSCS_CSE", "Te"), ("TIPI_Ag", "CFA"),
    ("TIPI_Ag", "Po"), ("TIPI_Co", "MRS_av"), ("TIPI_Ex", "Te"),
    ("TIPI_Ne", "MRS_co"), ("TIPI_Ne", "Po"), ("TIPI_Op", "MRS_in"),
    ("TIPI_Op", "SSCS_CPI"), ("TIPI_Op", "SSCS_CSE"), ("TIPI_Op", "TIPI_Ex"),
    ("Te", "TIPI_Co")
]

# NetworkXの有向グラフオブジェクトを作成
G = nx.DiGraph(directed_edges)

# サイクルが存在するかチェック (トポロジカルソートの前提条件)
if not nx.is_directed_acyclic_graph(G):
    print("エラー: グラフにサイクル（循環）が存在するため、階層を計算できません。")
    # サイクルを見つけて表示
    cycles = list(nx.simple_cycles(G))
    print("検出されたサイクル:", cycles)
else:
    print("--- 因果グラフの階層分析結果 ---")
    
    # Kahnのアルゴリズムを応用して階層を計算
    in_degree = {node: G.in_degree(node) for node in G.nodes()}
    
    level = 0
    levels = defaultdict(list)
    
    # 現在の階層のノードリストを初期化 (レベル0)
    current_level_nodes = [node for node, degree in in_degree.items() if degree == 0]
    
    processed_nodes_count = 0
    
    while current_level_nodes:
        # 現在のレベルのノードを結果に格納
        levels[level] = sorted(current_level_nodes)
        processed_nodes_count += len(current_level_nodes)
        
        next_level_nodes = []
        # 現在のレベルの各ノードについて処理
        for node in current_level_nodes:
            # そのノードから出ている先のノードのin-degreeを1減らす
            for neighbor in G.successors(node):
                in_degree[neighbor] -= 1
                # in-degreeが0になったら、次のレベルの候補に追加
                if in_degree[neighbor] == 0:
                    next_level_nodes.append(neighbor)
        
        # 次のレベルへ
        level += 1
        current_level_nodes = next_level_nodes

    # 結果の表示
    if processed_nodes_count != len(G.nodes()):
         print("警告: すべてのノードを処理できませんでした。グラフが連結していない可能性があります。")

    for lvl, nodes in sorted(levels.items()):
        print(f"\n[レベル {lvl}]")
        if lvl == 0:
            print("  (他の変数から影響を受けない、因果の開始地点)")
        
        # そのレベルのノードがどこに向かっているかを表示
        for node in nodes:
            targets = sorted(list(G.successors(node)))
            if targets:
                print(f"  - {node} --> {', '.join(targets)}")
            else:
                # どこにも向かっていないノードはゴール地点の候補
                print(f"  - {node} (この変数から出ていくパスなし)")

print("\n---------------------------------")
