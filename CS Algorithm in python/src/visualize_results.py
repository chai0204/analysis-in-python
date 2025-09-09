
import networkx as nx
import matplotlib.pyplot as plt
import japanize_matplotlib

# ユーザーから提供された分析結果データ
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

undirected_edges = [
    ("MW_S", "TIPI_Ne")
]

# グラフオブジェクトの作成
G = nx.DiGraph()
G.add_edges_from(directed_edges)
# 無向エッジのノードも追加（有向エッジに含まれていない場合のため）
for u, v in undirected_edges:
    G.add_node(u)
    G.add_node(v)

# グラフの可視化
plt.figure(figsize=(20, 16))

# ノードの配置を計算（バネモデル）
pos = nx.spring_layout(G, k=1.2, iterations=70, seed=42)

# ノードの描画
nx.draw_networkx_nodes(G, pos, node_size=3500, node_color='skyblue', alpha=0.9, linewidths=1, edgecolors='black')

# ラベルの描画
nx.draw_networkx_labels(G, pos, font_size=12, font_weight='bold')

# 有向エッジの描画
nx.draw_networkx_edges(G, pos, edgelist=directed_edges,
                       edge_color='#555555', width=1.5,
                       arrowstyle='-|>', arrowsize=25,
                       node_size=3500)

# 無向エッジの描画（赤色の破線）
nx.draw_networkx_edges(G, pos, edgelist=undirected_edges,
                       edge_color='red', width=2, style='dashed',
                       arrows=False,
                       node_size=3500)

# グラフのタイトルと設定
plt.title("因果分析結果の可視化グラフ", size=25)
plt.axis('off')
plt.tight_layout()

# ファイルに保存
output_filename = 'causal_graph.png'
plt.savefig(output_filename, format='png', dpi=300)

print(f"グラフが '{output_filename}' として保存されました。")
