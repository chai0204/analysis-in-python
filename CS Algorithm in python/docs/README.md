# CSアルゴリズムに基づく因果探索分析

## 概要

このリポジトリは、磯崎隆司氏の論文(2014)で提案されたCSアルゴリズム（Combining-Stage Algorithm）に基づき、観測データから変数間の因果構造を推定するためのPythonスクリプト群です。

分析の目的に応じて、2つの異なるアプローチを提供します。

1.  **無向グラフ分析 (`cs_algorithm_undirected.py`)**: 変数間の関連性の「骨格」を素早く把握したい場合に適しています。出力は変数間の関連性の有無とその強さを示します。
2.  **有向グラフ分析 (`cs_algorithm_directed.py`)**: 論文のロジックに完全に準拠し、因果の方向性を含んだ詳細なグラフ（PDAG）を構築します。

いずれのスクリプトも、分析プロセスを詳細にログ出力し、最終的な分析結果をJSONファイルとして`output`ディレクトリに保存します。

## ファイル構成

-   `src/`
    -   `cs_algorithm_undirected.py`: **無向グラフ分析**を実行するスクリプト。
    -   `cs_algorithm_directed.py`: **有向グラフ分析**を実行するスクリプト。
    -   `prepare_sachs_data.py`: Sachs(2005)のベンチマークデータを整形するための補助スクリプト。
-   `data/`
    -   `sample_data.csv`: 動作確認用のサンプルデータ（300件 x 20変数）。
-   `output/`
    -   分析結果のJSONファイルが保存されるディレクトリ。
-   `docs/`
    -   `README.md`: このファイル。

## 必要なライブラリ

-   pandas
-   pingouin
-   networkx

```bash
pip install pandas pingouin networkx
```

## データセットについて

### 1. サンプルデータ (`data/sample_data.csv`)

リポジトリに最初から含まれている、動作確認用のデータです。分析スクリプトがすぐに実行できるように、300件・20変数のランダムな整数データが格納されています。

### 2. Sachs (2005) ベンチマークデータ (要準備)

因果探索アルゴリズムの性能評価における標準的なベンチマークとして広く知られているデータセットです。細胞内のプロテインシグナル伝達に関するデータであり、専門家によって検証された「正解」の因果ネットワークが存在するため、アルゴリズムの妥当性を評価するのに適しています。

**注意**: このデータセットは著作権で保護されているため、リポジトリに直接含めることはできません。利用するには、以下の「実行方法」に従ってご自身で準備していただく必要があります。

## 実行方法

### ステップ1: データの準備

分析したいデータを選択・準備します。

-   **サンプルデータを利用する場合**: `data/sample_data.csv`が既にあるため、このステップは不要です。

-   **Sachsデータセットを利用する場合**:
    1.  **ダウンロード**: [こちらのURL](http://www.bnlearn.com/book-crc/code/sachs.interventional.txt.gz)から`sachs.interventional.txt.gz`をダウンロードします。
    2.  **展開と配置**: ダウンロードしたファイルを展開（解凍）し、中にある`sachs.interventional.txt`をこのプロジェクトの`data`ディレクトリに配置します。
    3.  **整形スクリプトの実行**: ターミナルで以下のコマンドを実行し、分析用の`data/sachs_data.csv`を生成します。
        ```bash
        python src/prepare_sachs_data.py
        ```

### ステップ2: 分析の実行

1.  **実行設定**: `src/cs_algorithm_undirected.py`または`src/cs_algorithm_directed.py`を開き、末尾の`main`関数内にある設定項目を編集します。

    ```python
    def main():
        # --- 設定項目 ---
        # 分析対象のCSVファイルパス
        INPUT_CSV_PATH = 'data/sample_data.csv'  # or 'data/sachs_data.csv'
        
        # 統計的検定の有意水準（α）
        SIGNIFICANCE_LEVEL = 0.05
        
        # 条件付き独立性検定で考慮する最大変数数
        MAX_CONTROL_VARS = 4
        
        # 分析結果を出力するJSONファイルのパス
        OUTPUT_JSON_PATH = 'output/causal_analysis_results_undirected.json'
    ```

2.  **スクリプト実行**: ターミナルで以下のいずれかのコマンドを実行します。

    ```bash
    # 無向グラフ分析の場合
    python src/cs_algorithm_undirected.py
    ```
    実行後、`output/causal_analysis_results_undirected.json` に結果が出力されます。

    ```bash
    # 有向グラフ分析の場合
    python src/cs_algorithm_directed.py
    ```
    実行後、`output/causal_analysis_results.json` に結果が出力されます。

### (補足) 他のスクリプトからの利用

各分析スクリプトは、関数として外部からインポートして利用することも可能です。

```python
from src.cs_algorithm_undirected import run_undirected_analysis
from src.cs_algorithm_directed import run_directed_analysis

# 無向グラフ分析の実行
run_undirected_analysis(
    input_csv_path='my_data.csv',
    significance_level=0.01,
    max_control_vars=5,
    output_json_path='output/my_undirected_results.json'
)

# 有向グラフ分析の実行
run_directed_analysis(
    input_csv_path='my_data.csv',
    significance_level=0.01,
    max_control_vars=5,
    output_json_path='output/my_directed_results.json'
)
```

## 参考文献

-   Isozaki, T. (2014). A Robust Causal Discovery Algorithm against Faithfulness Violation. *Information and Media Technologies*, 9(1), 121–131.
-   Sachs, K., Perez, O., Pe'er, D., Lauffenburger, D. A., & Nolan, G. P. (2005). Causal protein-signaling networks derived from multiparameter single-cell data. *Science*, 308(5721), 523-529.

---
# (English Version)

# Causal Discovery Analysis with the CS Algorithm

## Overview

This repository provides a set of Python scripts to infer causal structures from observational data, based on the Combining-Stage (CS) algorithm proposed by Isozaki (2014).

It offers two distinct approaches based on your analysis goals:

1.  **Undirected Graph Analysis (`cs_algorithm_undirected.py`)**: Suitable for quickly identifying the "skeleton" of relationships between variables. The output indicates the presence and strength of associations.
2.  **Directed Graph Analysis (`cs_algorithm_directed.py`)**: Fully adheres to the paper's logic to construct a detailed graph (PDAG) that includes causal directions.

Both scripts provide detailed logs of the analysis process and save the final results as a JSON file in the `output` directory.

## File Structure

-   `src/`
    -   `cs_algorithm_undirected.py`: Script for **undirected graph analysis**.
    -   `cs_algorithm_directed.py`: Script for **directed graph analysis**.
    -   `prepare_sachs_data.py`: A helper script to format the Sachs (2005) benchmark dataset.
-   `data/`
    -   `sample_data.csv`: A sample dataset for quick testing (300 samples x 20 vars).
-   `output/`
    -   The directory where the resulting JSON files are stored.
-   `docs/`
    -   `README.md`: This file.

## Requirements

-   pandas
-   pingouin
-   networkx

```bash
pip install pandas pingouin networkx
```

## About the Datasets

### 1. Sample Data (`data/sample_data.csv`)

This repository includes a sample dataset for immediate use, containing 300 samples of 20 random integer variables.

### 2. Sachs (2005) Benchmark Dataset (Requires Preparation)

This is a well-known benchmark dataset in causal discovery. As it comes with a ground-truth network, it is ideal for validating the algorithm's performance.

**Note**: Due to copyright, this dataset is not included directly. Please follow the steps below to prepare it.

## Usage

### Step 1: Prepare the Data

-   **To use the sample data**: No preparation is needed.

-   **To use the Sachs dataset**:
    1.  **Download**: Download the data from [this URL](http://www.bnlearn.com/book-crc/code/sachs.interventional.txt.gz).
    2.  **Extract & Place**: Unzip the file and place `sachs.interventional.txt` into the `data` directory.
    3.  **Run Prep Script**: Execute the following command to generate the formatted `data/sachs_data.csv`.
        ```bash
        python src/prepare_sachs_data.py
        ```

### Step 2: Run the Analysis

1.  **Configure**: Open `src/cs_algorithm_undirected.py` or `src/cs_algorithm_directed.py` and edit the configuration variables in the `main` function at the end of the script.

    ```python
    def main():
        # --- Configurations ---
        # Path to the input CSV file
        INPUT_CSV_PATH = 'data/sample_data.csv'  # or 'data/sachs_data.csv'
        
        # Significance level (alpha) for statistical tests
        SIGNIFICANCE_LEVEL = 0.05
        
        # Maximum number of control variables for conditional independence tests
        MAX_CONTROL_VARS = 4
        
        # Path for the output JSON file
        OUTPUT_JSON_PATH = 'output/causal_analysis_results_undirected.json'
    ```

2.  **Execute**: Run one of the following commands in your terminal.

    ```bash
    # For undirected graph analysis
    python src/cs_algorithm_undirected.py
    ```
    The results will be saved to `output/causal_analysis_results_undirected.json`.

    ```bash
    # For directed graph analysis
    python src/cs_algorithm_directed.py
    ```
    The results will be saved to `output/causal_analysis_results.json`.

### (Optional) Importing as a Module

You can also import and use the analysis functions in other scripts.

```python
from src.cs_algorithm_undirected import run_undirected_analysis
from src.cs_algorithm_directed import run_directed_analysis

# Run undirected analysis
run_undirected_analysis(
    input_csv_path='my_data.csv',
    significance_level=0.01,
    max_control_vars=5,
    output_json_path='output/my_undirected_results.json'
)

# Run directed analysis
run_directed_analysis(
    input_csv_path='my_data.csv',
    significance_level=0.01,
    max_control_vars=5,
    output_json_path='output/my_directed_results.json'
)
```

## References

-   Isozaki, T. (2014). A Robust Causal Discovery Algorithm against Faithfulness Violation. *Information and Media Technologies*, 9(1), 121–131.
-   Sachs, K., Perez, O., Pe'er, D., Lauffenburger, D. A., & Nolan, G. P. (2005). Causal protein-signaling networks derived from multiparameter single-cell data. *Science*, 308(5721), 523-529.
