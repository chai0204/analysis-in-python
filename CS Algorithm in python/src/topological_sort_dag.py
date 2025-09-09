
import networkx as nx
from collections import defaultdict
import json

def analyze_dag_levels(directed_edges: list, graph_name: str):
    """
    与えられた有向グラフの階層をトポロジカルソートで計算し、表示する関数。
    
    Args:
        directed_edges (list): (source, target) のタプルで構成される有向辺のリスト。
        graph_name (str): 分析対象を識別するための名前。
    """
    print(f"--- [{graph_name}] の階層分析結果 ---")
    
    G = nx.DiGraph(directed_edges)

    # 1. サイクル検出 (トポロジカルソートの前提条件)
    if not nx.is_directed_acyclic_graph(G):
        print("エラー: グラフにサイクル（循環）が存在するため、階層を計算できません。")
        try:
            cycles = list(nx.simple_cycles(G))
            print("検出されたサイクル (最初の5件まで):", cycles[:5])
        except Exception as e:
            print(f"サイクルの検出中にエラーが発生しました: {e}")
        print("-" * 30 + "\n")
        return

    # 2. 階層計算 (Kahnのアルゴリズム)
    in_degree = {node: G.in_degree(node) for node in G.nodes()}
    
    level = 0
    levels = defaultdict(list)
    
    # レベル0のノード (入次数が0のノード) から開始
    queue = [node for node, degree in in_degree.items() if degree == 0]
    
    processed_nodes_count = 0
    
    while queue:
        levels[level] = sorted(queue)
        processed_nodes_count += len(queue)
        
        next_queue = []
        for node in queue:
            for neighbor in sorted(list(G.successors(node))):
                in_degree[neighbor] -= 1
                if in_degree[neighbor] == 0:
                    next_queue.append(neighbor)
        
        queue = next_queue
        level += 1

    # 3. 結果表示
    if processed_nodes_count != len(G.nodes()):
         print("警告: すべてのノードを処理できませんでした。グラフが連結していない可能性があります。")

    for lvl, nodes in sorted(levels.items()):
        print(f"\n[レベル {lvl}]")
        if lvl == 0:
            print("  (因果の開始地点)")
        
        for node in nodes:
            targets = sorted(list(G.successors(node)))
            if targets:
                print(f"  - {node} --> {', '.join(targets)}")
            else:
                print(f"  - {node} (因果のゴール地点候補)")
    
    print("-" * 30 + "\n")


def main():
    # --- JSONファイルから有向パスを読み込む ---
    json_file_path = 'output/causal_analysis_results.json'
    directed_edges = []
    
    try:
        with open(json_file_path, 'r', encoding='utf-8') as f:
            results = json.load(f)
        
        for item in results:
            if item.get("向き") == "-->":
                directed_edges.append((item["変数1"], item["変数2"]))
        
        print(f"'{json_file_path}' から {len(directed_edges)} 件の有向パスを読み込みました。")

    except FileNotFoundError:
        print(f"エラー: JSONファイル '{json_file_path}' が見つかりません。")
        return
    except Exception as e:
        print(f"エラー: JSONファイルの読み込み中に問題が発生しました: {e}")
        return

    # --- ここまで ---

    # 読み込んだデータで実行
    analyze_dag_levels(directed_edges, "JSONからの分析結果")


if __name__ == '__main__':
    main()
