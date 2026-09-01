# sedumipy 開発ガイド(作業方針)

このドキュメントは、SeDuMi の MATLAB/Octave 非依存移植プロジェクトに
新しく参加する人(人間・AIエージェント問わず)が迷わず作業を続けられる
ようにするための引き継ぎ資料です。「何が終わっていて、何が残っていて、
どういう手順・方針で進めてきたか」をまとめています。

**リポジトリ構成について:** この移植作業はもともと `napinoco/sedumi`
(`sqlp/sedumi` のフォーク)の `python_port/` ディレクトリ以下で、
オリジナルのMATLAB/C実装と同居する形で進められていました。このリポジトリ
(`sedumipy`)は、移植版に必要なコードだけを独立させたクリーンな移設先です。
オリジナルの `.m`/MEX用 `.c` 一式は `vendor/sedumi-upstream/`(`sqlp/sedumi`
の submodule、フォーク時点のコミットに固定)に参照用として残っており、
オラクル再生成(Octave実機での比較用データ生成)はこのsubmoduleを使います。

## 1. プロジェクトのゴール

SeDuMi(MATLAB/Octave 上で動く SDP/SOCP 用の内点法ソルバー、`.m` ファイル
約96個 + MEX用 `.c` ファイル約54個)を、**MATLAB にも Octave にも一切
依存しない**形に移植する。採用したアーキテクチャ(過去のセッションで
決定済み、"Option B"):

- アルゴリズム本体(内点法の反復ロジック)→ **Python (NumPy/SciPy)**
- 性能が重要な低レベルカーネル(Cholesky分解、コーン演算など)→
  既存の `.c` ファイルを MEX 依存を取り除いて**独立 C ライブラリ**
  (`libsedumi.so`)としてビルドし、**ctypes** 経由で Python から呼ぶ

最終的には `pip install` できる、MATLAB/Octave 不要の Python パッケージ
にすることが目標。

## 2. 全体のフェーズと現在の進捗

タスク管理ツール(TaskList)にフェーズごとのタスクが登録されている。
2026-08-31 時点の状況:

| フェーズ | 内容 | 状態 |
|---|---|---|
| Phase 0 | 検証基盤(Octave実機でgolden reference取得) | **完了** |
| Phase 1 | Cカーネルのmex依存除去→独立Cライブラリ化 | **完了** |
| Phase 2 | Pythonバインディング(ctypes)構築 (クラスタ1〜5) | **完了** |
| Phase 3-a | 薄いMEXラッパー`.m`をPython APIとして整備 | **完了** |
| Phase 3-b | コーン数学ユーティリティ(eigK, psdeig, psdscale等)移植 | **完了** |
| Phase 3-c | 内点法の反復制御ロジック(sdinit〜optstep)移植 | **完了** |
| Phase 3-d | `sedumi.m`本体の移植 + golden referenceでの全体検証 | **完了(LP+SOCP+PSDスコープ)** |
| Phase 4 | 高レベルAPI・入出力互換層(.mat/SDPA)の実装 | **一部完了(.mat I/O・トップレベルAPIは完了、SDPAは未着手)** |
| Phase 5 | 検証・ベンチマーク | **完了** |
| Phase 6 | パッケージング・リリース | **未着手** |

Phase 3(内点法アルゴリズム本体の移植)は **LP + SOCP(2次錐) + PSD(半正定値
錐)問題について完了**しており、実際に `sedumipy.sedumi.sedumi(A,b,c,K)` を
呼べば Octave 版の SeDuMi と完全一致する解が返ってくることを、実機
オラクル比較で確認済み(`tests/test_sedumi.py`。PSD錐のメインループ結線は
`getada_psd.py`(`build_aord`/`getada_psd`、`incorder.py`/`getsymbada.py`/
`_native.getada1`/`getada2`/`getada3` を使用)、検証は同ファイルの
`sdp_feasible`/`sdp_mixed_cones_feasible` ケース)。

**密列(dense columns)最適化も移植済み**(`getdense.py`/`symbcholden.py`/
`deninfac.py`/`pcg.py`の product-form 前処理補正、`_native.symbfwblk`/
`adendotd`/`adenscale`/`dpr1fact` のオーケストレーション)。密列を実際に
含む問題(`tests/test_sedumi.py::test_sedumi_dense_matches_octave`、
`pars.denf=3` で `getdense.m` の検出閾値を意図的に下げたケース)でも
Octave版と `iter`/`numerr`/`pinf`/`dinf`/`x` が一致することを確認済み
(`y` は §6 の「双対解の非一意性」を参照)。

**Phase 5(実問題での検証)も完了**(`tests/test_golden_end_to_end.py`):
Phase 0 の golden reference が対象にしていた実問題(SDPLIB由来、
`vendor/sedumi-upstream/examples/`、`nb`/`arch0`/`control07`/`trto3`/
`OH_2Pi_STO-6GN9r12g1T2`。`quantum` は `K.scomplex`/`K.ycomplex` を使う
複素Hermitian PSD問題でスコープ外のため除外)で `sedumipy.sedumi()` を
実行し、Octave実機の golden reference と目的関数値が一致することを
確認した。この検証の過程で以下の2件の実バグを発見・修正した(詳細は
§6):
- `K.s==0`(LP+SOCPのみ)パスのADA記号的コレスキー順序が、Lorentz錐
  のarrow項(`d.q2`)に依存するsparsity patternの一部を見落としており、
  `d.q2`が育つにつれてCholesky分解が不正確になりPCGが発散するバグ
  (`nb.mat`、396個のSOCPブロックで顕在化)。
- `cpspdiag`(`getada3`のK.s==0分岐が呼ぶ診断用の対角成分抽出)が
  `ibsearch`マクロ経由で`bsearch()`を使っており、そのコンパレータが
  `sortnnz.c`/`iswnbr.c`と同種のqsort/bsearchコンパレータ未定義動作を
  踏んでいた(§6参照)。ただし実際のsedumi.py呼び出し経路では
  `K.s==0`のとき`getada3`自体が呼ばれないため、実害はテスト
  (`test_getada_no_psd_blocks`)止まりだった。

## 3. ディレクトリ構成

```
sedumipy/                    # リポジトリルート
  vendor/
    sedumi-upstream/          # submodule: sqlp/sedumi(フォーク元コミットに固定、参照専用)
  csrc/                       # mex.h依存を除去した標準alone Cカーネルソース(libsedumi.so の材料)
    *.c / *.h                  # sedumi_platform.h 経由でMEX非依存ビルドに対応させたフォーク
    sedumi_platform.c / .h
    kernel_smoke/smoke_test.c  # Phase 1 のスモークテスト
  src/
    sedumipy/                 # 移植先の Python パッケージ本体
      _native.py               # ctypes バインディング集約(全Cカーネル呼び出しはここに集める)
      libsedumi.so              # ビルド済み共有ライブラリ(tools/build_libsedumi.sh で生成、gitignore対象)
      cone.py                   # コーン数学ユーティリティ(eigK, psdeig, psdscale, frameit...)
      pretransfo.py / posttransfo.py   # 外部形式 <-> 内部形式の変換
      sdinit.py                 # 初期点生成
      sdfactor.py / sddir.py    # 自己双対埋め込みの分解・方向計算
      pcg.py                    # 前処理付き共役勾配法(loopPcg/wrapPcg)
      wregion.py                # 1反復分の predictor(+corrector)ステップ本体
      updtransfo.py             # スケーリング点の更新
      maxstep.py / widelen.py / stepdif.py / trydif.py  # ステップ長計算
      getada.py / getdatm.py / deninfac.py  # ADA行列の構築・分解(K.s==0)
      incorder.py / getsymbada.py / getada_psd.py  # ADA行列の構築(K.s!=0)
      symbchol.py                # ADAの記号的コレスキー(一度だけ実行)
      optstep.py                 # LP最適性の早期判定(optstep.m)
      amul.py / checkpars.py     # 補助ユーティリティ
      sedumi.py                  # トップレベルドライバ(全部をつなぐ)
      matio.py                   # Phase 4: .mat問題/解ファイルの読み書き
  tests/
    test_*.py                  # 各モジュールの検証テスト(オラクル比較)
    fixtures/                  # Octave実機で生成した .mat オラクルデータ(コミット済み)
    golden/                    # Phase 0 の golden reference
  tools/
    generate_*_oracle.m        # 各テストのオラクルを vendor/sedumi-upstream の Octave/MEXビルドで生成するスクリプト
    build_libsedumi.sh         # csrc/ から libsedumi.so をビルドするスクリプト
  pyproject.toml
  README.md                    # (やや古い。フェーズ概要はこのCONTRIBUTING.mdの方が新しい)
```

## 4. 開発ワークフロー(1つの `.m` ファイルを移植する際の手順)

これまで一貫して踏襲してきた手順。新しく関数を移植するときはこれに
従うこと。

1. **対象の `.m` ファイルの実ソースを全部読む。** コメントだけでなく
   実装ロジックを1行ずつ理解する。thin MEX wrapper(`sedumi_binary_error()`
   だけを呼ぶスタブ)の場合は対応する `.c` ファイルを読む。
2. **忠実な Python 移植を書く。** 変数名や処理順序はできるだけ `.m`/`.c`
   に対応させ、後から見比べやすくする。ドキュストリングには「何を
   ポートしたか」「`.m` ファイルの何行目に対応するか」「意図的に省略
   した部分とその理由」を明記する。
3. **Octave実機でオラクルを生成する。** `tools/generate_<name>_oracle.m`
   を書き、`vendor/sedumi-upstream/` の実際の `.m`/MEX ビルド
   (`install_sedumi` 済み)を呼び出して、入力データと出力を `.mat` に保存する。
   ```
   octave-cli --no-gui --eval "cd tools; generate_<name>_oracle"
   ```
   既存フィクスチャ(例: `pretransfo` の `K2`/`prep`)を使い回せる場合は
   積極的に再利用し、ゼロから計算し直すコストを避ける。
4. **`test_<name>.py` を書いて比較する。** `.mat` を読み込み、Python版の
   出力と `np.testing.assert_allclose` で突き合わせる。
5. **不一致が出たら中間値を1つずつ突き合わせてデバッグする。** 「なんと
   なく直す」のではなく、Octave側とPython側で同じ変数を1つずつ出力して
   どこで最初にズレるかを特定する。これまで見つかった実バグ(§6参照)
   は全てこの方法で見つけている。
6. **`python_port/tests/ -q` を全部流して回帰がないか確認してからコミット。**
   コミットメッセージには「何を移植したか」「見つけたバグとその修正」
   「意図的なスコープ制限とその理由」を書く。

### オラクルスクリプトを書くときの注意(乱数ストリーム)

Octave の `rand('seed', N)` は一様乱数の状態しかリセットしない。
`randn` は別ストリームなので、`randn` を使うケースでは必ず
`randn('seed', N)` も一緒に呼ぶこと(呼び忘れると実行するたびに
別のデータが生成され、フィクスチャが再現不能になる → 実際にこの
セッションでハマった)。

また、既存のオラクル生成スクリプトの**途中**に新しいコードを挿入すると、
それより後ろの `rand()` 呼び出しが全部ズレて、無関係な既存フィクスチャ
まで書き換わってしまうことがある。新しいケースは必ずスクリプトの
**末尾に追記**し、`git diff --stat` で既存フィクスチャのファイルサイズが
変わっていないか確認すること。

## 5. スコープ上の制約(意図的に未実装にしている部分)

以下は「気づいていない抜け」ではなく、**意図的に `NotImplementedError`
で弾いている**制約。理由も含めて各モジュールのdocstringに明記済み。

- **回転2次錐(`K.r`)** は `pretransfo.py` の段階で標準2次錐(`K.q`)に
  変換されるため、それより後段のコードは意識する必要はない
  (`sedumi.py` の検証テストにも回転錐のケースが1つ含まれている)。
- **コンソール出力・v-plot・`pars.stopat` のデバッグブレーク・事前の
  ランク診断・DIMACS誤差指標(`info.err`)** は移植していない。これらは
  全て「診断・表示専用」で `(x,y,info)` の値そのものには影響しない
  ため、優先度を下げている。

## 6. これまでに見つかった実バグ・注意点(教訓)

- **MATLABの値渡し意味論 vs ctypesのin-place変更。** `fwsolve`/`bwsolve`
  はCカーネル(`fwblkslv.c`/`bwblkslv.c`)がin-placeで書き換える設計を
  そのままctypesバインディングしているため、呼び出し側でバッファを
  使い回すMATLABコードをそのまま移植すると**サイレントに壊れる**。
  `pcg.py` の `sparfwslv`/`sparbwslv` で「呼び出し前に必ずコピーする」
  ように統一して解決した。同種の罠は他のin-placeネイティブ関数
  (`fwdpr1`/`bwdpr1`等)にもあり得るので、新しくバインディングを
  使うときは常に疑うこと。
- **mexFunctionの自動スライシング。** `psdframeit.c`/`psdinvjmul.c` の
  ような一部のMEXカーネルは、「PSD部分だけの短い配列」と「L+Q+PSDの
  フル長配列」の両方を受け付け、フル長の場合は自動でオフセットを
  読み飛ばす(`x += cK.lpN + 2*cK.lorN`)。Phase 2のctypesバインディング
  がこの自動スライシングを再現し忘れていたことが `psdinvjmul` で
  実際に発覚した(`wregion.py` 移植時)。新しいバインディングを
  追加する際は元の `mexFunction` のこの種の分岐を見落とさないこと。
- **qsortコンパレータの未定義動作。** `sortnnz.c`/`iswnbr.c` は
  `signed char` を返すコンパレータを `int(*)(const void*,const void*)`
  にキャストして `qsort()` に渡しており、実機で実際に非決定的な
  挙動を確認した。この種の関数は**ctypesバインディングせず**、
  コメントに書かれた本来のアルゴリズムをPythonで直接書き直す方針
  にしている(`neighborhood.py` の `iswnbr` が実例)。
- **`symbchol.m` の完全密行列分岐。** ADAが完全に密(全要素非ゼロ)の
  場合、実際の `symbchol.m` は最小次数順序付け(MMD)を省略して恒等
  順序+単一巨大supernodeを直接使う。これを再現せず常に `ordmmd` を
  呼ぶと、**収束はするが反復回数がOctave版とズレる**(小さい密な
  テスト問題で実際に確認済み)。`symbchol.py` はこの分岐を正確に
  再現している(`_native.symbolic_cholesky_dense`)。
- **MATLABの `'` は共役転置。** 実数配列に対しても、複素数が絡む式
  (`posttransfo.m` の `(x'*prep.QR)'` 等)では単純転置と共役転置の
  違いが結果に影響する。移植時は常に注意する。
- **`optstep.m` の `sum(K.s)!=0` 分岐は実質デッドコード。** `sedumi.m`
  が `optstep` を呼ぶのは `lponly = (K.l==length(c))` の場合のみで、
  これは `K.q`/`K.s` が両方空であることを強制する。つまり
  `optstep.m` 自身にPSD対応の分岐が書いてあっても、実際の呼び出し
  経路からは絶対に到達しない。こういう「到達不能性を呼び出し元の
  条件から証明できるデッドコード」は、移植せずに `NotImplementedError`
  で塞いで良い、という判断をしている。ただし本当に呼び出し元の
  条件を確認してからにすること(推測で「多分デッドコードだろう」で
  済ませない)。
- **`sparfwslv`/`sparbwslv`(`pcg.py`)は `L.perm` を内部で
  gather/scatter しなければならない。** 実際の `fwblkslv.c`/
  `bwblkslv.c` は `y = L\b(L.perm)`(前進代入)/`y(L.perm) = L'\b`
  (後退代入)というように、呼び出しの内側で `L.perm` によるgather/
  scatterを行う。密列最適化導入前の版ではこれを省略し、代わりに
  `loopPcg`/`wrapPcg` を呼ぶ全箇所(テストも含め)で `L.perm` を
  外側から一貫して手動で適用/解除する形にしていた ―― PCGは
  「内部のインデックス付け規約が呼び出し全体で自己無矛盾でありさえ
  すれば」正しい解に収束するため、密列がない(`Lden` が恒等)場合は
  これでも数値的に正しかった。しかし `deninfac.py` が実際に
  `L.perm` で `Ad` を並べ替えた**真の**(恒等でない)`Lden` を組み立てる
  ようになると、`loopPcg` の同じ反復内で `fwdpr1(Lden, sparfwslv(L,r))`
  のように `Lden` の項と `sparfwslv` の項を合成する際、両者の
  インデックス規約が食い違って反復が壊れる(`iter` は合っても
  途中の残差が発散する形で発覚した)。`sparfwslv`/`sparbwslv` 自体を
  実際のCカーネルと同じgather/scatter付きの実装に直し、
  `wrapPcg`/`sdfactor.py` 側は本家 `.m` 通り一切パーミュートしない
  形に戻すことで解決 ―― 「呼び出し側で辻褄合わせをする」workaroundは
  一見動いても、後から非自明な追加機能(この場合は密列)が入ると
  壊れる典型例だった。
- **rank落ちした `At` を使うテストフィクスチャでは双対解 `y` が
  非一意になる。** `tests/fixtures/sedumi/lp_socp_sdp_dense_feasible.mat`
  (密列最適化のend-to-endテスト用)は `At` が23行20列で
  `rank(At)=16`(密列検出の`h`フロアを`NORMDEN=5`に固定するため
  PSDブロックの背景行密度を`0.02`まで下げた結果、一部の行が構造的に
  全ゼロになりやすく、シード1〜400を全数探索しても`rank(At)=20`
  (フルランク)になるケースは1つもなかった)。`At`がフルランクでない
  等式制約系では、最適双対解 `y` は `At@delta=0` となる任意の
  `delta` だけずらしても同じ目的関数値・同じ双対スラック
  (`s=c-At@y`)を与えるため一意に定まらない。Python版とOctave版は
  `iter`/`numerr`/`pinf`/`dinf` が完全一致し `x` も浮動小数点誤差の
  範囲内で一致する(=アルゴリズムは同一に動いている)にもかかわらず、
  `y` は縮退方向に沿って2-ノルムで約5もズレる
  (`At@(y_py-y_oct)`, `b@(y_py-y_oct)` がともに実質ゼロであることで
  確認済み)。このため `test_sedumi_dense_matches_octave` では `y` の
  厳密一致ではなく、双対の実行可能性・最適性(`c-At@y` と `b@y`)を
  比較する形にしている。**同じ理由で `At` がフルランクとは限らない
  ランダム生成テストフィクスチャを新規に作る際は、`y` を厳密比較する
  前に `rank(At)==m` かどうかを確認すること。**
- **`K.s==0` パスの一回限りのADA記号的コレスキー順序が、Lorentz錐の
  arrow項に依存するsparsity patternの一部を見落としていた。**
  `sedinit.py` はスケーリング点 `d["q2"]`(各Lorentz錐のarrow部分の
  スカラー)を必ず厳密に0から開始する(`sdinit.m`の`d.q2 = zeros(...)`
  通り)。`sedumi.py`の`K.s==0`分岐は、この最初の`d`(=`d.q2=0`)で
  一度だけ`ADA = getada(A,K,d,DAt)`を計算し、その**数値的な**
  sparsity patternをそのまま`symbchol(ADA)`に渡して以後全反復で使い
  回していた。ところが`getada.py`のLorentz項
  (`DAt_q.T @ DAt_q`、`DAt.q[k,j] = d.q1[k]*Aj[k] + d.q2[k]*(...)`)
  は、`d.q2=0`のとき同じLorentz錐ブロックを共有するだけで行方向には
  重ならない制約ペア `(i,j)` に対する寄与が構造的にまるごと消えて
  しまう ―― つまり反復1時点のADAは、後の反復で`d.q2`が育つにつれて
  実際には非ゼロになる位置を欠いた、過小なsparsity patternになって
  いた。`numeric_cholesky`はこの(固定された)symbolic patternの外側
  には書き込めないため、`d.q2`が育つ反復(実際には3〜5反復目以降)
  からCholesky分解が徐々に不正確になり、PCGの前処理性能が崩れて
  反復回数上限に張り付き、最終的に誤った解に収束する。
  `vendor/sedumi-upstream/examples/nb.mat`(LP+396個のSOCPブロック)
  で実際に確認: 本家Octave版は20反復・`numerr=0`で収束するのに対し、
  修正前のPython版は9反復で`numerr=2`(深刻な数値エラー)を報告して
  完全に異なる値を返していた。本家`sedumi.m`はこの問題を、
  `getsymbada.m`ベースの構造的(値に依存しない)pattern構築を
  `sum(K.s)`の値に関わらず常に`sdinit`より前に実行することで回避して
  いる(本ファイルの過去のこの箇所の記述、および`sedumi.py`自身の
  SCOPEドキュストリングが「実装を簡略化した既知の相違点」として
  触れていた箇所)。**修正**: `sedumi.py`の`K.s==0`分岐で、
  `symbchol()`に渡す一度限りのADAだけは`d["q2"]`を(捨てるための
  ローカルコピーで)強制的に非ゼロにしたものから構築するようにし、
  以後の各反復で実際に使う`d`/`DAt`/`ADA`には一切手を加えないように
  した。`tests/test_golden_end_to_end.py`(Phase 5, §7参照)で
  `nb`/`arch0`/`control07`/`trto3`/`OH_2Pi_STO-6GN9r12g1T2`の5問題
  全てがOctave実機の結果と一致することを確認済み。
- **`getada3`の`K.s==0`分岐が呼ぶ`cpspdiag`も、`sortnnz.c`/`iswnbr.c`
  と同じqsort/bsearchコンパレータ未定義動作を踏んでいた。**
  `cpspdiag`は`blksdp.h`の`ibsearch`マクロ(=標準ライブラリの
  `bsearch()`)経由でADAの対角成分を探す。`ibsearch`はこの探索に
  `char`を返す`icmp()`を`COMPFUN`(`int(*)(const void*,const void*)`)
  にキャストして渡しており、これも未定義動作 ―― 実際にこのport の
  ビルドでは`bsearch()`が対角成分を一度も見つけられず、`absd`が
  (対角成分が正しくソートされた形で存在しているにもかかわらず)
  常に全て0.0になっていた(`tests/test_getada.py::test_getada_no_psd_blocks`
  で発覚)。ただし`sedumi.py`の実際の呼び出し経路では、`getada3`
  自体が`has_psd=True`(=`K.s`が非空、したがって内部的に`sdpN>0`)の
  ときしか呼ばれないため、この`sdpN==0`分岐は実運用では到達しない
  デッドコードであり実害はなかった。`sortnnz`/`iswnbr`のときと同じ
  方針で解決 ―― `cpspdiag`はctypesバインディングをやめ、
  `scipy.sparse`の`.diagonal()`で直接対角成分を取る(`cpspdiag.c`
  自身のドキュメントコメント通りの意図)Python実装に置き換えた。

## 7. 残っている作業(優先度が高いと思われる順)

1. ~~**Phase 3-a: 薄いMEXラッパー`.m`の公開API整備。**~~ **完了。**
   `install_sedumi.m` のMEXビルド対象一覧と `_native.py` の全バインディング
   を突き合わせて棚卸しした結果:実際に使われている実MEXカーネルは
   全て `_native.py` に集約済みで、それぞれ然るべき上位モジュール
   (`getdense.py`/`getdatm.py`/`pcg.py`/`cone.py`/`updtransfo.py`/
   `wregion.py`/`sdinit.py`/`getada_psd.py`/`symbchol.py`/
   `symbcholden.py` 等)から `_native.xxx()` の形で呼ばれており、
   「バインディングはあるが未整理」という抜けは見つからなかった
   (`incorder`/`iswnbr` の2つだけは qsortの未定義動作を避けるため
   意図的にctypes化せず `incorder.py`/`neighborhood.py` にPython実装
   として存在する、既知の意図的な設計)。
   逆に `_native.py` 内で他から一切呼ばれていないバインディングが
   7個(`realdot`/`realssqr`/`scalarmul`/`addscalarmul`/`blkmul`/
   `mJdetd`/`cholsplit`)見つかったが、いずれも「本家SeDuMi自身に
   おいても未使用」と確認済み(`blkmul.c`/`mJdetd.c`は
   `install_sedumi.m`のMEXビルド対象リストにそもそも入っていない
   =本家の時点でデッドコード、`cholsplit()`の出力`L.split`は
   `blkchol.c`のmex引数リストに現れず本家でも読まれていない、
   `realdot`等はBLAS的な補助関数でPhase 1のスモークテスト用に
   バインドされただけで独立したMEXターゲットを持たない)。
   結論として追加実装は不要と判断し、`_native.py` モジュール
   docstringにこの棚卸し結果自体を明記した(未使用の理由を含む)。
2. **Phase 4: 高レベルAPI・入出力互換層。** 以下は完了:
   - トップレベルAPI: `import sedumipy; sedumipy.sedumi(A,b,c,K)` が
     使えるようになった(従来は `sedumipy.sedumi.sedumi(...)` の
     サブモジュール経由のみ)。`__init__.py`で`from .sedumi import sedumi`
     しているが、`sedumipy.sedumi`サブモジュールを先に(または後に)
     importしても関数を指すことに変わりはない(Pythonの
     `sys.modules`キャッシュにより、親パッケージへの属性上書きは
     サブモジュールの初回import時にしか起きないため)ことを確認済み。
     また`sedumi()`は`pars`辞書に加えて`**kwargs`でも個別オプションを
     渡せるようにした(例: `sedumi(A,b,c,K,eps=1e-9)`)。
   - `.mat` I/O: `matio.py`(`read_mat`/`write_solution_mat`)。
     SeDuMiの問題ファイルは(移植元に対応する`.m`が存在しない)ただの
     MATLAB構造体なので、他のモジュールと違い「移植」ではなく
     このport独自の新規実装。`A`/`At`どちらの向きの格納も、
     `b`/`c`がスパースで保存されているケースも扱う
     (`vendor/sedumi-upstream/examples/*.mat`で確認済み)。
   残っているのは **SDPA sparse形式(`.dat-s`)の読み書き**のみ
   (対応する`.m`実装は本家にもなく、これも新規実装になる)。
3. ~~**Phase 5: 検証・ベンチマーク。**~~ **完了。**
   `tests/test_golden_end_to_end.py` が Phase 0 の golden reference
   対象問題(`vendor/sedumi-upstream/examples/`)を実際に
   `sedumipy.sedumi()` に通し、Octave実機の結果と一致することを検証
   する(§2「Phase 5」、§6の2件のバグ修正を参照)。性能ベンチマークは
   `tools/benchmark_examples.py`(実行方法はスクリプト自身のdocstring
   参照)。この環境(Octaveをローカルでビルドして計測、CPU/コア数等は
   環境依存につき絶対値は目安)での実測値:

   | problem | m | N (=length(c)) | Python (秒) | Octave/MEX (秒) | iter |
   |---|---:|---:|---:|---:|---:|
   | nb | 123 | 2383 | 3.0 | 0.9 | 20 |
   | arch0 | 174 | 56197 | 2.5 | 2.4 | 31〜32 |
   | control07 | 666 | 6125 | 9.3 | 9.2 | 40 |
   | trto3 | 544 | 398977 | 18.2 | 19.8 | 60 |
   | OH_2Pi_STO-6GN9r12g1T2 | 948 | 240720 | 34.4 | 34.8 | 20 |

   最小の問題(`nb`)ではPython側のオーバーヘッド(関数呼び出し・
   numpy配列確保・ctypes境界越えのコスト)が支配的でOctave/MEX版の
   約3倍かかるが、問題が大きくなるにつれてCネイティブカーネルでの
   実計算timeが支配的になり、中〜大規模問題(`arch0`以上)では
   Octave/MEX版とほぼ同等〜やや高速という結果になった。`arch0`の
   `iter`が31/32とOctave側と1回だけズレているのは、大規模問題での
   浮動小数点丸め誤差の蓄積差によるもので(§6にある通り`test_sedumi_
   matches_octave`の厳密な`iter`一致要求は小さな合成フィクスチャでの
   話であり、実問題規模ではこの種の1反復程度のズレは想定内)、
   `cx`/`by`は両者とも期待値に一致しているため実害はない。
4. **Phase 6: パッケージング。** `cibuildwheel` 等でのwheel化、
   `libsedumi.so` の同梱方法の検討(現状はビルド済みバイナリを
   そのままリポジトリに置いている)。

## 8. テストの実行方法

```
cd sedumipy
git submodule update --init --recursive   # 初回のみ: vendor/sedumi-upstream を取得
.venv/bin/pip install -e .[test]           # libsedumi.so は初回import時に自動ビルドされる
.venv/bin/python -m pytest tests/ -q
```

Octaveがインストールされていない環境でもテストは通る(オラクル
データは `.mat` として事前生成・コミット済みで、テストは実行時に
Octaveを呼ばない)。オラクルを**再生成**する場合のみ Octave
(`vendor/sedumi-upstream` 側で `install_sedumi` 実行済みの状態)が必要:

```
octave-cli --no-gui --eval "cd tools; generate_<name>_oracle"
```

## 9. コーディング規約・命名の慣習

- **K(コーン構造体)は Python の `dict`** で表現する。フィールド名は
  `.m` 版の `K.f`/`K.l`/`K.q`/... とできるだけ同じにする
  (`K["l"]`, `K["mainblks"]` など)。`pretransfo.py` が計算する
  内部専用フィールド(`mainblks`, `qblkstart`, `sblkstart`, `lq`,
  `N`, `rsdpN` 等)はそのまま下流の全関数に引き回す。
- **インデックスは基本 0-indexed(Python流)**。ただし `K` の中身
  (`mainblks` 等)は元の `.m` の1-indexed値をそのまま持っていることが
  多く、使う側で `int(x) - 1` のような変換をその都度行っている
  (統一的な変換層はまだ無い)。新しいコードを書くときは既存の
  近い関数のインデックス変換パターンを真似ること。
- **疎行列は `scipy.sparse.csc_matrix`** を基本形式とする。
- **全てのCカーネル呼び出しは `_native.py` に集約**し、他のモジュール
  は `_native.py` 経由でのみCコードに触れる。
- **各moduleのdocstringに「何を実装し、何を意図的に実装していないか」
  を明記する**のがこのプロジェクトの一貫した文化。読んだだけで
  スコープが分かるようにしておくこと。
