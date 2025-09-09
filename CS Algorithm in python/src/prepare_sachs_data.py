# prepare_sachs_data.py
import re
import os

# --- 設定 ---
# 入力ファイル名（bnlearn.comからダウンロードした元のファイル）
INPUT_FILENAME = "sachs.interventional.txt"
# 出力ファイル名（分析スクリプトで使用するCSVファイル）
OUTPUT_FILENAME = "sachs_data.csv"

def clean_value(value):
    """ハイフン付きの値を処理する。例: '1-1' -> '1'"""
    return value.split('-')[0]

def main():
    """スクリプトのメイン処理"""
    print(f"--- {INPUT_FILENAME} の処理を開始します ---")

    # 入力ファイルの存在チェック
    if not os.path.exists(INPUT_FILENAME):
        print(f"\nエラー: '{INPUT_FILENAME}' が見つかりません。")
        print("Sachs et al. (2005)のデータセットを以下のURLからダウンロードし、")
        print("このスクリプトと同じディレクトリに配置してください。")
        print("URL: http://www.bnlearn.com/book-crc/code/sachs.interventional.txt.gz")
        print("(ダウンロード後、.gzファイルを展開して .txt ファイルを配置してください)")
        return

    try:
        with open(INPUT_FILENAME, 'r', encoding='utf-8') as f_in:
            lines = f_in.readlines()

        # ヘッダー行の処理: 引用符を削除し、スペースをカンマに置換
        header = lines[0].strip().replace('" ", ', '').replace('"', '')

        # データ行の処理
        cleaned_lines = [header]
        for line in lines[1:]:
            if not line.strip(): continue # 空行はスキップ
            # 1つ以上のスペースを区切り文字として値を分割
            values = re.split(r'\s+', line.strip())
            # 各値をクリーンアップ
            cleaned_values = [clean_value(v) for v in values]
            cleaned_lines.append(",".join(cleaned_values))

        # CSVファイルとして書き出し
        with open(OUTPUT_FILENAME, 'w', encoding='utf-8') as f_out:
            f_out.write("\n".join(cleaned_lines))
        
        print(f"\n処理が完了しました。")
        print(f"クリーンなデータが '{OUTPUT_FILENAME}' として保存されました。")

    except Exception as e:
        print(f"\nエラーが発生しました: {e}")

if __name__ == '__main__':
    main()
