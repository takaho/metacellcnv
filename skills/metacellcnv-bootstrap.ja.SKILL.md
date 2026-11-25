---
name: metacellcnv-bootstrap
description: >
  CellRangerの出力ディレクトリだけを渡されたとき、metacellcnvパイプライン
  (https://github.com/takaho/metacellcnv)をGitHubから取得し、適合するPython環境と
  依存パッケージを一式インストールし、ゲノム名からGTF/染色体マップを自動取得し、
  前処理→metacell構築→CNV推定→悪性判定→DE→可視化レポートまでを実行して結果を返す。
  ユーザー側に事前のソフトウェアインストールは一切要求しない。CellRanger出力から
  CNV解析結果が欲しい、metacellcnvを動かしたい、腫瘍細胞をCNVで同定したい、
  といった依頼で起動する。
---

# metacellcnv ブートストラップ & 実行スキル

あなた(エージェント)は、ユーザーがCellRangerの出力しか持っていない環境で、
`metacellcnv` パイプライン一式を自力で用意し、最後まで実行する。ユーザーに
「まず○○をインストールしてください」とは絶対に言わない。インストールが
必要なものは全部あなたがやる。

## 0. 受け取るパラメータ

ユーザーの依頼文から以下を読み取る。明示が無いものは「既定値」を使い、
実行前に「この値を使う」と一行で伝える(聞き返して手を止めない)。

| パラメータ | 必須 | 説明 | 既定値 |
|---|---|---|---|
| `cellranger_dir` | ○ | `filtered_feature_bc_matrix/` を含むCellRanger出力ディレクトリ(複数可、サンプルごとに1つ) | — |
| `sample_id` | — | 各 `cellranger_dir` に対応するサンプル名 | ディレクトリ名から推定 |
| `genome` | ○(`gtf`未指定時) | 参照ゲノム名。例: `GRCh38`, `GRCm39`, `CanFam3.1`, `CanFam6`, `Dog10K_Boxer_Tasha`、または任意のNCBIアセンブリ名/アクセッション | — |
| `gtf` | — | 既存のGTFファイルパス。指定時はゲノム名からの自動取得をスキップ | 無し(自動取得) |
| `chromosome_map` | — | NCBI assembly report / 2列の染色体変換表 | GTF取得元に応じて自動決定 |
| `mito_chromosome` | — | mtDNAの配列ID | GTFから自動検出 + 既知ID |
| `markers` | — | 細胞型マーカー。`dog`/`human`/`mouse` またはCSVパス | `genome`から推定(推定不能ならユーザーに確認) |
| `cnv_method` | — | `infercnv`(既定, 窓移動平均)/ `bins` | `infercnv` |
| `cnv_reference` | — | `celltype`(腫瘍組織+正常細胞混在) / `none`(純粋な細胞株) / `external` | サンプル構成から推定(判断材料が無ければ`celltype`) |
| `sample_kind` | — | `tissue` / `cell-line` | `auto`(パイプライン既定) |
| `n_cpus` | — | 並列数 | マシンの論理コア数 |
| `out_dir` | — | 成果物の出力先 | `./metacellcnv_run/<timestamp>/` |
| `repo_url` | — | metacellcnvのGitHub URL | `https://github.com/takaho/metacellcnv.git` |
| `repo_ref` | — | ブランチ/タグ/コミット | リポジトリ既定ブランチ(通常 `main`) |
| `python_version` | — | 使用するPythonバージョン | `3.11`(リポジトリ `INSTALL.md` 推奨。`>=3.10`なら可) |

## 1. 環境診断(ユーザーの既存環境には触れない)

1. `python3 --version` に加え、`pyenv`, `uv`, `conda`/`mamba`, `asdf` の有無を
   `command -v` で確認する。
2. 要求Pythonバージョン(既定3.11、`>=3.10`で可)を用意できる手段を選ぶ優先順位:
   `uv`(`uv venv --python 3.11`)→ `conda`/`mamba`(`conda create -n metacellcnv python=3.11`)
   → `pyenv`(無ければそのツール自体をユーザーの許可を得てインストール)→
   システムの`python3`が既に`>=3.10`ならそれで`venv`を作る。
3. いずれにせよ **専用の仮想環境**(名前例: `metacellcnv`)を作り、システムの
   Python/パッケージは一切変更しない。既に同名の環境があれば使い回す(再構築しない)。
4. 仮想環境のパスと選んだ方法を一行で報告する。

## 2. リポジトリ取得

1. 作業用ディレクトリ(例: `./_metacellcnv_src/`、ユーザーの `cellranger_dir` や
   出力先とは分離する)を用意する。
2. 既にそこにクローン済みなら `git fetch && git checkout <repo_ref>` で更新、
   無ければ `git clone --depth 1 -b <repo_ref> <repo_url> .` で取得する
   (ブランチ未指定ならdepth 1の通常clone)。
3. `git`が使えない/ネットワークがgit向けに閉じている場合は、GitHubのtarball API
   (`https://github.com/<owner>/<repo>/archive/refs/heads/<ref>.tar.gz`)を取得して展開する。
4. 取得できたら `README.md` と `INSTALL.md` を読み、CLIオプションが本スキル記載の
   ものと食い違っていないか(バージョンが進んでいる可能性がある)軽く確認する。

## 3. 依存パッケージのインストールと自己検証

1. 仮想環境内で、リポジトリの `INSTALL.md` に従いインストールする
   (最小構成ではなく「本番相当」一式を入れる):
   ```
   pip install numpy pandas scipy scikit-learn anndata scanpy \
               pydeseq2 scrublet plotly
   ```
   可視化を省略したい等の指示が無ければ `plotly` も入れる。
   `--cnv-method bins` や旧matplotlibレポートを明示的に求められた場合のみ
   追加パッケージ(`infercnvpy`, `matplotlib`)を足す。
2. **自己検証を必ず行う**: `python smoke_test.py` を実行し、
   `ALL SMOKE TESTS PASSED` が出ることを確認してから先に進む。
   失敗したら本番データには触れず、エラー内容をそのままユーザーに報告して止まる
   (憶測で回避策を入れない)。
3. この時点でインストールと自己検証が完了した一行サマリをユーザーに出す
   (例: 「環境構築OK(uv, python 3.11, metacellcnv HEAD, smoke test PASSED)」)。

## 4. ゲノムアノテーション(GTF)の自動取得

`gtf` が明示されていればここは全てスキップしてそれを使う。

1. **キャッシュ確認**: `reference/<genome>/` のような専用キャッシュディレクトリを
   作業域に用意し、既に同じ `genome` 名で取得済みのGTF(+assembly report)が
   あれば再ダウンロードしない(2回目以降の実行を速くする)。
2. **ヒト/マウスはGENCODE優先**(染色体名が `chr1...` で、`--chromosome-map` 不要になる):
   - `genome` が `GRCh38`/`hg38`/`human` 系 →
     `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/`
     のディレクトリ一覧を取得し、`gencode.v*.annotation.gtf.gz` を取得。
   - `genome` が `GRCm39`/`mm39`/`mouse` 系 →
     `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/latest_release/`
     から同様に取得。
   - 旧アセンブリ(`GRCh37`/`hg19`, `GRCm38`/`mm10`)を明示された場合は
     GENCODEの該当アーカイブ版ディレクトリを探す。
   - 取得後、`gtf-gene-id` は `gene_name`(= `--gtf-gene-id auto` で自動検出されるはず。
     念のため `--preflight-only` で確認)。
3. **それ以外(イヌを含むヒト・マウス以外、または特定アクセッション指定)は
   NCBI RefSeq/GenBankアセンブリ経由**:
   - `genome` が既にNCBIアクセッション(`GCF_*`/`GCA_*`)ならそのまま使う。
   - アセンブリ名(例 `CanFam6`, `Dog10K_Boxer_Tasha`, `CanFam3.1`)の場合は
     NCBI Datasets API / `datasets` CLI、または Entrez `esearch -db assembly`
     で最新アクセッションを解決する。**正確な対応表をこちらで決め打ちしない**
     — 必ずその場でAPI照会して確認する(アセンブリは改版されるため)。
     このプロジェクトでは過去に `CanFam3.1` = `GCF_000002285.5` を使った実績がある
     (参考情報であり、最新アセンブリが欲しい場合はAPIの結果を優先する)。
   - 解決したアクセッションから、NCBI genomesのFTP/HTTPSパス規則
     (`.../genomes/all/GCF/NNN/NNN/NNN/GCF_XXXXXXXXX.X_<AssemblyName>/`)
     の下にある `*_genomic.gtf.gz` と `*_assembly_report.txt` を取得する
     (`datasets download genome accession <ACC> --include gtf,seq-report` が使えるなら
     それで一括取得してよい)。
   - 取得した `*_assembly_report.txt` を `--chromosome-map` としてそのまま渡す
     (RefSeqの `NC_...` を `chr1...` に変換するため)。この場合 `--gtf-gene-id gene`
     になる(GENCODEの `gene_name` とは違う属性名)。
4. ミトコンドリア配列IDは、既知の対応(イヌCanFam系列は `NC_002008.4` など)が
   使えればそれを、無ければGTFからの自動検出(パイプライン既定動作)に任せる。
5. 取得物のパスと由来(URL・アクセッション)を出力ディレクトリの
   マニフェストに書き残す(再現性のため)。

## 5. 入力検証(preflight)

各 `cellranger_dir` について、本実行の前に必ず一度:

```
python metacellcnv.py --cellranger-dir <dir> --gtf <gtf> --gtf-gene-id auto \
    [--chromosome-map <map>] --markers <markers> --preflight-only
```

を実行し、遺伝子名とGTFの対応率・mtDNA遺伝子の有無・マーカー遺伝子の有無・
染色体名を確認する。対応率が低い、マーカーがほぼ当たらない、といった警告が
出た場合は、`--gtf-gene-id` の値やゲノム名の取り違えを疑い、原因を特定してから
本実行に進む(警告を無視して続行しない)。

## 6. 本実行(サンプルごとに3段階)

```
# 6a. 前処理
python metacellcnv.py --scanpy \
    --cellranger-dir <cellranger_dir> --sample-id <sample_id> \
    --gtf <gtf> --gtf-gene-id auto [--chromosome-map <map>] \
    --mito-chromosome <mito> \
    --out-dir <out_dir>/prep/<sample_id>

# 6b. metacell -> CNV -> 悪性判定 -> DE
python metacellcnv.py \
    --cellranger-dir <cellranger_dir> --sample-id <sample_id> \
    --gtf <gtf> --gtf-gene-id auto [--chromosome-map <map>] \
    --mito-chromosome <mito> \
    --markers <markers> \
    --cnv-method <cnv_method> --cnv-reference <cnv_reference> \
    [--cnv-external-reference <path> ] [--sample-kind <sample_kind>] \
    --n-cpus <n_cpus> \
    --out-dir <out_dir>/results/<sample_id>

# 6c. 可視化レポート
python metacellcnv.py --visualize \
    --results-dir <out_dir>/results/<sample_id> \
    --prep-dir <out_dir>/prep/<sample_id>
```

複数サンプルがあり、かつ同種の正常組織サンプルを外部参照に使いたい場合
(`--cnv-reference external`)は、`build_external_cnv_reference.py` 相当の
プロファイル化ステップを先に済ませてからこの手順に入る。

## 7. 進捗報告とETA

各段階の開始時に一行で「今どの段階か」を報告する(例:
`[2/3 samples] Case1: metacell構築中(seacells-fw-iters 50) 開始 14:32`)。
ETAは次の方針で出す — **最初から確度の高い数字を装わない**:

1. 各段階の開始・終了時刻を自分で記録する。
2. 最初に終わる軽い段階(前処理 `--scanpy`)の「経過時間 / 細胞数」を
   その場のマシンでの実測スループットの基準値として使う。
3. 重い段階(metacell構築・CNV推定)は、同じ細胞数・遺伝子数に対して
   この基準値から線形外挿してETAを出し、「前段の実測値からの概算」である旨を明記する。
4. 標準出力に途中経過(反復回数、処理済みバーコード数など)が出るスクリプトは
   それを読み、一定間隔(例: 30秒〜1分、または出力の区切りごと)でETAを更新する。
   既存の出力に進捗情報が全く無い段階は、「進捗不明。経過時間のみ表示」と正直に言う。
5. サンプルが複数ある場合は、既に完了したサンプルの所要時間を次サンプルのETAに使う
   (1サンプル目が最も不確実、2サンプル目以降は精度が上がる)。
6. 想定より大幅に(目安2倍以上)時間がかかっている段階は、固まっていないか
   (CPU使用率・ログの更新)を確認し、一言添えてから報告を続ける。

## 8. 失敗時の扱い

- `smoke_test.py` が失敗 → 本番データに進まず、エラーをそのまま報告して停止。
- `--preflight-only` で対応率が極端に低い/マーカーがほぼ無い →
  `--gtf-gene-id` やゲノム取り違えを疑い、原因を特定してから再実行。
  原因不明のまま本実行に進まない。
- ネットワークが閉じていてGTF取得やリポジトリ取得ができない →
  その旨と、必要なURL/アクセッションを明示してユーザーに手動配置を依頼する
  (無言で諦めない)。

## 9. 完了後

- サンプルごとに `prep/<sample_id>/`, `results/<sample_id>/`, 可視化レポート
  (自己完結HTML)のパスを一覧で示す。
- 使った実際のパラメータ(ゲノム・アクセッション・GTF由来・各CLIオプション)を
  `out_dir` 直下にマニフェストとして残し、再実行時に差分だけ分かるようにする。
- 仮想環境・クローン済みリポジトリ・ダウンロード済みゲノムはそのまま残し、
  次回以降の呼び出し(別の `cellranger_dir` だけが変わるケースなど)では
  環境構築・リポジトリ取得・GTF取得を全てスキップして本実行だけ行う
  (= 2回目以降は大幅に速くなることをユーザーに伝えてよい)。
