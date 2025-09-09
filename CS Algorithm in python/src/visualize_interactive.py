import networkx as nx
from pyvis.network import Network
import json

# --- JSONファイルから分析結果を読み込む ---
json_file_path = 'output/causal_analysis_results.json'
directed_edges_data = []
undirected_edges_data = []

try:
    with open(json_file_path, 'r', encoding='utf-8') as f:
        results = json.load(f)
    
    for item in results:
        u = item.get("変数1")
        v = item.get("変数2")
        direction = item.get("向き")
        weight = item.get("偏相関係数")

        if not all([u, v, direction, weight is not None]):
            print(f"警告: 不完全なデータ項目をスキップしました: {item}")
            continue

        if direction == "-->":
            directed_edges_data.append((u, v, weight))
        elif direction == "---" or direction == "<-->":
            undirected_edges_data.append((u, v, weight))

    print(f"'{json_file_path}' からデータを正常に読み込みました。")
    print(f"  - 有向パス: {len(directed_edges_data)}件")
    print(f"  - 無向・双方向パス: {len(undirected_edges_data)}件")

except FileNotFoundError:
    print(f"エラー: JSONファイル '{json_file_path}' が見つかりません。")
    exit()
except json.JSONDecodeError:
    print(f"エラー: JSONファイル '{json_file_path}' の形式が正しくありません。")
    exit()
except Exception as e:
    print(f"エラー: データの読み込み中に問題が発生しました: {e}")
    exit()


# --- Pyvisによるグラフ生成 ---

# Pyvis Networkオブジェクトの初期化
net = Network(height="800px", width="100%", bgcolor="#222222", font_color="white", notebook=True, directed=True)

# ノードとエッジをグラフに追加
all_nodes = set()
for u, v, _ in directed_edges_data:
    all_nodes.add(u)
    all_nodes.add(v)
for u, v, _ in undirected_edges_data:
    all_nodes.add(u)
    all_nodes.add(v)

for node in all_nodes:
    net.add_node(node, label=node, size=25, font={'size': 20})

for u, v, weight in directed_edges_data:
    color = 'red' if weight < 0 else 'cyan'
    net.add_edge(u, v, value=abs(weight), title=f"Strength: {weight}", color=color)

for u, v, weight in undirected_edges_data:
    color = 'red' if weight < 0 else 'magenta'
    net.add_edge(u, v, value=abs(weight), title=f"Strength: {weight}", dashes=True, color=color)

# 物理演算のON/OFFを切り替えるUIコントロールを追加
net.show_buttons(filter_=['physics'])

# HTMLファイルとして出力（新しいファイル名）
output_filename = 'output/interactive_causal_graph_from_json.html'
net.show(output_filename)

print(f"\nインタラクティブグラフが '{output_filename}' として保存されました。")
print("このHTMLファイルをWebブラウザで開いてください。")
print("\n【操作方法】")
print("グラフの下部に表示される操作パネルから 'physics' -> 'enabled' のチェックを外すと、ノードの動きが止まり、自由に配置できます。")
