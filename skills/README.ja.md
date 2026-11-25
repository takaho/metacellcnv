# metacellcnv-bootstrap スキルの使い方

このディレクトリの `metacellcnv-bootstrap.ja.SKILL.md` / `metacellcnv-bootstrap.en.SKILL.md` は、
CellRanger の出力ディレクトリだけを渡せば、環境構築(Python venv 作成・依存パッケージ導入)
から [metacellcnv](https://github.com/takaho/metacellcnv) の取得、ゲノム名からのGTF/染色体マップ
自動取得、パイプライン本実行、進捗・ETA報告までをLLMエージェントに一気通貫でやらせるための
「スキル(エージェント向け指示書)」です。日本語版と英語版は内容は同一です。

> **検証状況について**: このスキルは **Claude(Claude Code / Cowork)でのみ動作確認済み** です。
> 下記の Codex / Antigravity / OpenCode での利用方法は、各ツールの一般的な「永続指示」の
> 仕組みに基づく案内であり、**このスキル自体をそれらのツールで実際に動かして確認したもの
> ではありません**。各ツールの挙動は変化が速いため、使う前に必ず各ツールの公式ドキュメント
> も確認してください。

## Claude での使い方(動作確認済み)

**Claude Code**:
```
~/.claude/skills/metacellcnv-bootstrap/SKILL.md        # 個人用(全プロジェクト共通)
<project>/.claude/skills/metacellcnv-bootstrap/SKILL.md # プロジェクト専用
```
のいずれかに、`metacellcnv-bootstrap.ja.SKILL.md` または `.en.SKILL.md` を
`SKILL.md` という名前でコピーして置いてください(日本語・英語を両方使いたい場合は
別ディレクトリ名、例 `metacellcnv-bootstrap-ja/` `metacellcnv-bootstrap-en/` に分けてください)。
配置後は「CellRangerの出力からCNVを出して」のように自然に話しかければ自動的に起動します。

**Cowork(このチャット)**: ダウンロードしたファイルをそのまま会話に貼り付けてもよいですし、
組織の設定によっては、ファイルカードから保存する/アップロードすることでスキルとして
登録できる場合があります。

## Codex CLI(OpenAI)での使い方 (未検証)

Codex は、リポジトリ直下に置く `AGENTS.md`(常時読み込まれる永続的なプロジェクト指示)と、
`~/.codex/prompts/` (ユーザー共通) または `.codex/prompts/` (プロジェクト専用) に置く
Markdown ファイルをスラッシュコマンド化する「カスタムプロンプト」機能を持っています。

- **常時有効にしたい場合**: `metacellcnv-bootstrap.en.SKILL.md`(英語版を推奨。Codexは
  英語指示の方が安定する傾向があります)の本文(冒頭の `---` で囲まれた YAML 部分は
  削除してよい)を `AGENTS.md` に貼り付ける、または追記する。
- **必要な時だけ呼び出したい場合**: 本文を `~/.codex/prompts/metacellcnv-bootstrap.md`
  (または `.codex/prompts/metacellcnv-bootstrap.md`)として保存し、`/metacellcnv-bootstrap`
  のように呼び出す。

正確なディレクトリ名・仕様は変更される可能性があるため、`codex --help` や公式ドキュメント
(developers.openai.com/codex)で最新情報を確認してください。

## Google Antigravity での使い方 (未検証)

Antigravity にはプロジェクト/ワークスペース/ユーザー単位の「Rules」という永続指示の仕組みが
あります。本文をそのまま Rules ファイルとして登録するか、不明な場合は新しいエージェント
セッションの最初のメッセージとして `metacellcnv-bootstrap.en.SKILL.md` の本文をそのまま
貼り付けてタスクを依頼してください。Rules の正確な置き場所は変更が比較的早い領域なので、
Antigravity 公式ドキュメント(antigravity.google/docs/rules)を確認してください。

## OpenCode での使い方 (未検証)

OpenCode も `AGENTS.md`(プロジェクト直下 / グローバル)による永続指示に対応しているほか、
`.opencode/command/*.md` に置くとスラッシュコマンド化される「Commands」、
`.opencode/agent/*.md` に置く「Agents」という、Markdown + フロントマター形式の仕組みを
持っています。

- **常時有効**: 本文を `AGENTS.md` に追記。
- **呼び出し式**: `.opencode/command/metacellcnv-bootstrap.md` として保存し
  `/metacellcnv-bootstrap` で呼び出す。

## どのツールでも使える最終手段

上記のいずれの仕組みも使えない/分からない場合でも、`metacellcnv-bootstrap.en.SKILL.md`
(または `.ja.`)の本文は普通の自然言語の指示書なので、シェルとファイル操作ができる
エージェントであれば、**新しいタスクの最初のメッセージとしてファイルの中身をそのまま
貼り付けて依頼する**だけで同じ手順を実行させられるはずです。冒頭の YAML フロントマター
(`---` で囲まれた `name`/`description`)はClaude固有の記法なので、無視されても問題ありません。
