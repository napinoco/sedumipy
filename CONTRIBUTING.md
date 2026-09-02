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
| Phase 4 | 高レベルAPI・入出力互換層(.mat/SDPA)の実装 | **完了** |
| Phase 5 | 検証・ベンチマーク | **完了** |
| Phase 6 | パッケージング・リリース | **一部完了(wheelビルド自体は動作確認済み、manylinux/CI上での検証は未着手)** |

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
      sdpa.py                    # Phase 4: SDPA sparse(.dat-s)形式の読み書き
  tests/
    test_*.py                  # 各モジュールの検証テスト(オラクル比較)
    fixtures/                  # Octave実機で生成した .mat オラクルデータ(コミット済み)
    golden/                    # Phase 0 の golden reference
  tools/
    generate_*_oracle.m        # 各テストのオラクルを vendor/sedumi-upstream の Octave/MEXビルドで生成するスクリプト
    build_libsedumi.sh         # csrc/ から libsedumi.so をビルドするスクリプト
  pyproject.toml
  setup.py                      # Phase 6: wheelビルド時にlibsedumi.soをコンパイルするbuild_extフック
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

- **本家の「全ブロック一括(all-or-nothing)分岐」は移植せず、ブロック
  ごとにクランプする。** `widelen.m`/`trydif.m`/`maxstep.m` の3箇所に、
  Lorentz錐ブロックの判別式に対する

  ```matlab
  tmp = halfxz.^2 - detxz;
  if all(tmp > 0)          % ← 全ブロックまとめて1回だけ判定
      lab2q = halfxz + sqrt(tmp);
  else
      lab2q = halfxz;      % ← 1ブロック巻き添えで全ブロックが劣化
  end
  ```

  という**グローバルな全か無かの分岐**が存在する(`maxstep.m` は
  `norm2` に対する同型のもの)。一見「`sqrt`に負値を渡さないための
  安全策」に見えるが、判別式は**厳密算術では必ず非負**である:
  この式が生む2つの固有値は `lab2q` と `detxz/lab2q` で、積が `detxz`、
  和が `2*halfxz` なので `tmp` は恒等的に `((lab1-lab2)/2)^2`、すなわち
  完全平方。負になるのは1ブロックの2固有値がほぼ一致したときの丸め
  誤差だけであり、しかも**完全に一致すると厳密に0**になって、本家の
  strictな `> 0` はこれも弾く(1つのLorentzブロックの複製で構成される
  構造的な問題では日常的に起きる)。実測でも、発火時の値は微小な負数
  ではなく**厳密に `-0.0`** だった。

  したがって移植側は**ブロックごとにクランプ**する:

  ```python
  lab2q = halfxz + np.sqrt(np.maximum(tmp, 0.0))
  ```

  判別式が非正なそのブロックについては本家のフォールバックと厳密に
  同値(`sqrt(0)==0`)、他の全ブロックについては本家が捨てていた正確な
  式を保ち、`sqrt`に負値を渡さないという本家の安全性もそのまま満たす
  ―― つまりトレードオフではなく**本家のどちらの分岐よりも厳密に正確**。

  `maxstep.m` のものは精度ではなく**安全性のバグ**である点に注意:
  フォールバック時に `norm2` は平方根を取らないまま `reltr - norm2` に
  使われる(二乗量を線形量から引く次元的な不整合)。判別式が1未満
  ―― スケーリング後は普通のこと ―― では `v < sqrt(v)` なので、
  本家版は錐の境界までのステップ長を**過大評価**する。nb_L2で計測
  すると64回中5回発火し、5回とも過大評価だった。

  効果(DIMACS、3箇所とも本家挙動 vs クランプで切り替えて計測):

  | 問題 | 本家の分岐 | クランプ | 実機Octave/MEX |
  |---|---|---|---|
  | nb_L2 | numerr=2, iter=10 | **numerr=0, iter=16** | numerr=0, iter=16 |
  | nql180old | numerr=2, iter=12 (cx=18.08 vs by=7.08) | **numerr=1, iter=42** (cx≈by 8桁) | numerr=1, iter=54 |
  | qssp30old | numerr=2 (cx=6.6017 vs by=6.3582) | **numerr=1** (cx=6.496695, 公表値6.4966749) | numerr=2(本家も失敗) |

  教訓として: **本家の分岐条件が「数学的にありえない場合」に対する
  防御である場合、それを逐語移植すると防御が過剰に広く効いてしまう
  ことがある。** 「本家にそう書いてあるから」は移植の既定方針として
  正しいが、その条件が守ろうとしている不変量(ここでは「判別式は
  完全平方だから非負」)を確認すれば、本家より厳密に良い実装が
  一意に決まることがある。
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
2. ~~**Phase 4: 高レベルAPI・入出力互換層。**~~ **完了。**
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
   - SDPA sparse形式(`.dat-s`)の読み書き: `sdpa.py`(`read_sdpa`/
     `write_sdpa`)。`read_sdpa`は`conversion/fromsdpa.m`の忠実な移植
     (Octave実機オラクルと一致確認済み、
     `tools/generate_sdpa_oracle.m`/`tests/fixtures/sdpa/`)。
     `write_sdpa`は本家に対応物がない新規実装(本家の
     `conversion/writesdp.m`はSDPA形式ではなく別形式のSDPpackを書き出す
     もので無関係)だが、実際に`vendor/sedumi-upstream/examples/
     arch0.mat`を`write_sdpa`で書き出し、それを実機Octaveの
     `fromsdpa.m`で読み戻して元の`(At,b,c)`と完全一致することを
     手動で確認済み(K.q/K.rはSDPA形式で表現できないため
     `write_sdpa`は明示的に`ValueError`で拒否する)。
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
4. **Phase 6: パッケージング。** ~~`libsedumi.so`の同梱方法の検討
   (現状はビルド済みバイナリをそのままリポジトリに置いている)~~ ――
   この記述自体が古かった: 実際には`libsedumi.so`/`.dylib`は
   `.gitignore`対象で**リポジトリにはコミットされておらず**、
   `_native.py`の`_ensure_built()`が初回import時に
   `tools/build_libsedumi.sh`を自動的に呼んでその場でビルドする
   (開発時の`pip install -e .[test]`はこれに依存している)、という
   のが実態だった。この方式はeditable installでは動くが、
   本物のwheelとしてインストールされた場合には壊れる
   (`csrc/`/`tools/`は`sedumipy`パッケージの一部として同梱されて
   いないため、importのたびにコンパイルし直すこともできないし、
   そもそもエンドユーザーの環境にgcc/BLASの開発ヘッダーが入っている
   保証もない)。

   **今回やったこと:** `setup.py`に`build_ext`をオーバーライドする
   カスタムステップ(`BuildLibsedumi`)を追加し、`pip install`/
   `python -m build --wheel`時に`tools/build_libsedumi.sh`と同じ
   コンパイルコマンドを1回だけ実行して`libsedumi.so`をビルド、
   `build_lib/sedumipy/`配下に直接配置することでwheelに同梱される
   ようにした(`_native.py`の`_ensure_built()`はそのままなので、
   wheelから入れた場合は既にファイルがあるため何もせず、editable
   installの場合は従来通り初回import時にビルドする、という二重の
   動作を両立させている)。ソースを持たない`Extension`を1つ
   登録しているのは、setuptoolsに「このwheelはプラットフォーム
   依存(`py3-none-any`ではない)」と正しく認識させるためだけの
   トリック。

   **この環境で実際に確認したこと:** `python -m build --wheel`で
   `cp311-cp311-linux_x86_64`タグ付きのwheelがビルドされ
   `libsedumi.so`が同梱されていること、そのwheelを(このリポジトリの
   `csrc`/`tools`に一切アクセスできない)独立した仮想環境に
   `pip install`し、`import sedumipy; sedumipy.sedumi(...)`が
   正しく動作することを確認した(このLinux環境限定)。

   **まだ検証できていないこと(この環境にDockerデーモンが無く
   実行できなかった):** `cibuildwheel`自体の実行(`pyproject.toml`
   に`[tool.cibuildwheel]`の設定は追加したが、manylinuxコンテナ上での
   実際のビルドは未検証)、macOS/Windowsでのビルド(`tools/
   build_libsedumi.sh`はgcc前提でWindowsのcl.exeには未対応)、
   `libblas`への動的リンクによる配布可搬性の問題(`ldd`で確認した
   限り`libblas.so.3`/`libopenblas.so.0`に動的リンクされており、
   本当にPyPI配布可能なmanylinux wheelにするには`auditwheel repair`
   でこれらを同梱するか静的リンクに切り替える必要がある。今回は
   その作業までは行っていない)。
5. **`getdatm.py`のOOM修正(`DAt.q`常時sparse化)が、逆にdenseな方が
   速い小〜中規模問題を遅化させていた件。修正済み。** `has_psd=False`
   (LP+SOCPのみ、K.s==0)経路で`DAt.q`/`ADA`を常にsparseで組み立てる
   ように直した結果(OOM対策としては正しい)、mが小さくADAが
   実質denseになる問題(例: `nb.mat`, m=123)では逆に約36%遅化して
   いた(sparse-sparse積`csr_matmat`が全体の56%を占めることを
   `cProfile`で確認)。`getsymbada()`が一度だけ計算する構造的ADA
   パターンの密度(既存の0.9閾値をそのまま流用)から`is_dense`を
   一度だけ判定し、`getDAtm()`/`getada()`がそれに応じてdense
   (numpy, BLAS matmul)/sparse(scipy, OOM回避)を切り替えるように
   `sedumi.py`から配線した。両分岐は値としてbug-for-bug同一(既存
   テスト・ベンチマーク全て回帰なしを確認済み)。`nb.mat`で実測:
   sparse固定(修正直後)2.12秒→ハイブリッド化後1.75秒(修正前の
   dense固定1.56秒に近い水準まで回復)。大規模問題(`nql180`/
   `qssp180`, m~1.3e5)は引き続きsparse経路を通るためOOMは再発しない。
6. ~~**DIMACS `nb_L2`のnumerr=2、根本原因を特定(修正は見送り)。**~~
   **解決済み。修正した(§6の「全ブロック一括分岐」の項を参照)。**
   前回セッションは原因を`widelen.py`の`all(tmp>0)`グローバル分岐まで
   完全に特定した上で「アルゴリズム自体のカオス的感度であり、分岐の
   numericsを変えるのは影響が読めない」として意図的に見送っていたが、
   **この判断は保守的すぎた**。`tmp`は任意の量ではなく`((lab1-lab2)/2)^2`
   という**完全平方**なので厳密算術では必ず非負であり、負になるのは
   丸め誤差だけ ―― つまり「どちらの分岐を選ぶべきか」は本来一意に
   決まる。しかも今回計測し直すと、問題の値は前回想定されていた
   `±1.78e-15`のような微小な負数ではなく**厳密に`-0.0`**で、本家の
   strictな`> 0`がゼロを弾いていただけだった(1つのLorentzブロックの
   複製で構成される問題では2固有値が厳密に一致するのは日常的)。
   `np.sqrt(np.maximum(tmp, 0.0))`とブロックごとにクランプすれば、
   当該ブロックには本家のフォールバックと厳密に同値、他の838ブロック
   には本家が捨てていた正確な式が残り、`sqrt`に負値を渡さない安全性も
   保たれる。結果:**nb_L2は numerr=2/iter=10 から numerr=0/iter=16 に
   改善**し、実機Octave/MEXビルドの反復回数(16)と公表値
   (`-1.62897198`、この移植版は`-1.628971959`)の両方に一致した。
   `stepdif=1`を強制する対症療法も不要になった(pars既定値は変更して
   いない)。同じ分岐が`trydif.m`(逐語コピー)と`maxstep.m`(同型、
   ただしフォールバックがより悪質)にもあったので3箇所とも修正済み。
   以下は原因特定に至るまでの前回セッションの記録(そのまま残す):

   実機Octave/MEXビルド(`vendor/sedumi-upstream`、この環境に
   `octave`/`liboctave-dev`/`libopenblas-dev`を追加導入して
   `install_sedumi -rebuild`でビルド)と、この移植版の両方から
   反復1〜3時点の`ADA`/`d`/`DAt.q`を書き出して直接突き合わせた結果:
   `d.l`/`d.det`(LP・trace部分)は反復3まで浮動小数点誤差レベル
   (~1e-13)で一致、`d.q1`/`d.q2`(Lorentz錐のスケーリング点)も
   反復2開始時点までは同様に一致するが、反復2のステップが生成する
   `d`(=反復3で使われる`d`)で`d.q1`の最大成分が約15%相対誤差で
   食い違う。これがちょうど`err["kcg"]`/`Lsd["kcg"]`が実機の1/1から
   6/5へ跳ね上がる反復と一致する ―― ここまでは前回セッションの記録。

   **今回、この続きを追って原因を完全に特定した。** 前回時点で
   「`updtransfo.py`は`updtransfo.m`と一行ずつ突き合わせ済みで差分
   なし」と判定されていたが、今回は行レベルの目視監査に加えて
   **実際に動かして検証**した: 実機Octave側の`wregion.m`内に一時
   デバッグ用`save()`を挿入し(コミットしない一時パッチ)、反復2の
   `wregion`が返す`xscl`/`zscl`/`w`(および直前の`d`/`K`)をそのまま
   `.mat`にダンプ、Pythonの`updtransfo()`にそれを**そのまま**渡して
   実行したところ、実機の反復3の`d.q1`/`d.q2`とビット単位で完全一致
   した(`max|q1|`小数第10桁まで完全一致)。つまり`updtransfo.py`は
   本当に無罪 ―― 差分は`updtransfo`より手前、`xscl`/`zscl`/`w`
   自体の計算にあることが確定した。

   次に、この移植版が自分自身で計算した反復2の`xscl`/`zscl`/`w`を
   同様にダンプして実機の値と直接比較したところ:`xscl`/`zscl`は
   絶対誤差~1.5e-13(浮動小数点ノイズレベル、4196次元ベクトルに対して
   2つの独立実装が一致する限界としては極めて良好)、`w["tdetx"]`/
   `w["tdetz"]`も同様に~1e-13〜1e-12で一致するにもかかわらず、
   **`w["lab"]`だけが最大絶対誤差7.6という大きさで食い違っていた**。
   `w["lab"]`は`tdetx`/`tdetz`からほぼそのまま計算される量のはずなので、
   この不釣り合いが決定的な手がかりになった。

   `widelen.py`の`_build_w()`(`widelen.m`の該当部分をそのまま移植した
   箇所)を見ると、Lorentz錐の固有値相当量`lab2q`は

   ```python
   tmp = halfxz**2 - detxz
   if np.all(tmp > 0):        # widelen.m: if all(tmp > 0)
       lab2q = halfxz + np.sqrt(tmp)
   else:
       lab2q = halfxz          # 839ブロック全部がこのフォールバックになる
   ```

   という、**839個のLorentz錐ブロック全体に対する単一のグローバルな
   all-or-nothing分岐**になっている(1ブロックでも判別式`tmp`が
   非正なら、他の838ブロックも含めて全部が精度の低いフォールバック式
   になる)。実際に反復2の`xscl`/`zscl`から`tmp`を計算してみると、
   839要素中838番目のブロック(0-indexで396番目)の`tmp`が
   実機ビルドでは`+1.78e-15`、この移植版では`-1.78e-15`と、
   **符号だけが反転するギリギリの値**になっていた(残り838ブロックの
   `tmp`はどちらも十分正)。`xscl`/`zscl`自体は~1e-13レベルで一致して
   いるのに、この1ブロックだけがちょうどゼロをまたぐ位置にあった
   ため、2つの独立した浮動小数点パイプライン(NumPy/SciPy+自前Cカー
   ネル vs Octave+実機BLAS、内部の総和順序やBLAS実装が異なる)の
   ごく僅かな丸め誤差の違いだけで`all(tmp>0)`の真偽が分かれ、
   その結果`lab2q`(ひいては`w["lab"]`全体)が全く別の式で計算される
   ことになり、以降の`updtransfo`でのスケーリング点更新が大きく
   分岐してしまう ―― というのが今回突き止めた完全な因果連鎖。
   (実際に実機の`xscl`/`zscl`をこの移植版の`_build_w()`にそのまま
   与えると`w["lab"]`はビット単位で実機と完全一致し、逆にこの移植版
   自身の`xscl`/`zscl`を与えると`tmp[396]`が負に転じて`all(tmp>0)`が
   Falseになることも確認済み。)

   **この`all(tmp>0)`というグローバル分岐自体は`widelen.m`本家に
   そのまま存在する設計**(1ブロックでも判別式が負になり得るなら、
   `sqrt`にNaN/複素数を渡すリスクを避けるため全ブロックを一括で
   保守的なフォールバック式に倒す、という意図的な安全策と読める)
   であり、この移植が独自に導入した誤りではない。2つの独立実装が
   ちょうどゼロをまたぐ量について反対側の丸め誤差を持つ、という
   状況そのものは(適切な乱数シードで再現できる)本質的にアルゴリズム
   側のカオス的感度であり、`updtransfo.py`/`widelen.py`/`tdet`/`ddot`
   等のどの一行を直しても解消しない類のものと判断する。したがって
   **意図的に修正を見送る**(§7旧項目6で既に触れていた「`stepdif=1`を
   強制すればnb_L2は解けるが、それはpars既定値を全問題に対して変更
   する副作用が大きいので採用しない」という判断と同じ理由 ―― 症状に
   対する場当たり的な修正ではなく、原因を完全に特定した上で「今の
   実装のままで良い」と判断した点が今回のセッションの進捗)。
7. ~~**`nql180`/`qssp180`のnumerr=2、再検証したら実は直っていた
   (`nql180old`は実質的なロバスト性のギャップとして依然未解決、
   `qssp180old`は検証完了・移植バグではないと確認)。**~~
   **`nql180old`のギャップも解決済み**(この項目末尾の「追記」を参照。
   §6の「全ブロック一括分岐」修正による)。以下は経緯の記録:
   前回セッションの記録(「OOM修正後も数反復でnumerr=2」)は古い情報
   だった。実際に上記5.のdense/sparseハイブリッド化後の版で
   `matio.read_mat()`経由で読み込んで実行したところ:
   - `nql180`(m=226,802、PSDブロックなし): **numerr=0, iter=16, 約39秒**
     で正常収束(DIMACS READMEの参照値が"N/A"のため、`cx`≈`by`・
     `feasratio`→1・`r0`=1e-8という内部無矛盾性で確認)。
   - `qssp180`: **numerr=0, iter=42, 約249秒**で正常収束(同様に
     参照値なしのため内部無矛盾性で確認)。
   一方、`nql30old`/`qssp30old`と同じ「old(旧式・非推奨)formulation」
   系列の`nql180old`は、この環境で実機Octave/MEXビルド
   (`install_sedumi -rebuild`、§6同様)と突き合わせたところ、
   **どちらも同じようには失敗しない**:実機は`iter=54`まで粘って
   `numerr=1`(精度は`pars.bigeps`止まりだが破綻はしない、コンソールに
   `skip=5361`という大量のCholeskyピボットスキップが出るほど数値的に
   厳しい問題ではある)で終えるのに対し、この移植版は`iter=27`
   (`feasratio=0.90`, `r0=0.53`)で`numerr=2`(完全失敗)を返す ――
   `nql30old`/`qssp30old`のような「本家でも同じく失敗するので移植バグ
   ではない」という単純な話ではなく、本家より早く・悪く失敗している
   という**実質的なロバスト性のギャップ**が残っている。`qssp180old`
   (このファミリーで最大、~36MB)は前回セッションの調査時間予算
   (Python版・実機版とも550秒)内にどちらも完走せず未検証だったが、
   **今回のセッションでより大きな時間予算を与えて両方とも完走させ、
   決着させた**: 実機Octave/MEXビルドは**iter=30、numerr=2**で失敗
   (`install_sedumi -rebuild`済みの環境で`tic`/`toc`実測1705秒)。
   この移植版も同じファイルを`matio.read_mat()`経由で読み込んで実行
   したところ**iter=30、numerr=2**で失敗(実測3557秒、反復ごとの
   経過時間を`wregion()`にモンキーパッチして記録し、順調に反復が
   進み続けていて途中でハングしていないことも確認済み)。**失敗する
   反復番号(30)まで完全に一致**しており、`nql180old`のような
   「本家より早く・悪く失敗する」ロバスト性のギャップは見られない
   ―― `qssp30old`/`nql30old`と同じ、本家でも同じように失敗する
   ジャンルの問題であることが確定した(移植バグではない)。
   **追記:`nql180old`のロバスト性ギャップも解消した(§6・上記項目6の
   `all(tmp>0)`修正による)。** この項目が「本家より早く・悪く失敗する
   =実質的なロバスト性のギャップ」として残していた唯一の未解決事項
   だったが、3箇所の分岐をブロック単位クランプに直した状態で計測し
   直すと:

   | nql180old | numerr | iter | cx vs by |
   |---|---|---|---|
   | 本家の分岐(3箇所とも) | **2**(完全失敗) | 12 | 18.08 vs 7.08 |
   | ブロック単位クランプ | **1** | 42 | 8桁一致 |
   | 実機Octave/MEXビルド | 1 | 54 | ― |

   本家挙動では反復12で `cx`/`by` が2.5倍も乖離したまま破綻するのに
   対し、クランプ版は `cx=0.9311428505`/`by=0.9311428684` と8桁一致
   まで収束して`numerr=1`で終える ―― **実機ビルド(iter=54)より
   少ない42反復で同じ`numerr=1`に到達**しており、「本家より早く・悪く
   失敗する」ギャップは解消した(DIMACS READMEの参照値は"N/A"のため
   内部無矛盾性で確認)。この問題ではフォールバックの発火率が異常に
   高く(`widelen`は12回中6回=50%、`maxstep`は50回中5回)、本家の
   分岐が常時暴発していたことが分かる。
   `qssp30old`も同修正で`numerr=2`→`numerr=1`になり、しかも
   **DIMACS READMEの公表値`6.4966749`に一致する解**(`cx=6.496695`)を
   返すようになった ―― この問題は**実機Octave/MEXビルド自身が
   `numerr=2`で失敗する**ので、ここは「本家と同等」ではなく
   **本家超え**である。`nql30old`(公表値`0.9460`)と併せて、
   `tests/test_benchmarks.py`の除外リストから公表値付きの
   パラメトライズドテストへ昇格させた。
8. **`qssp180old`のcProfileから発覚したPythonレベルの性能バグ2件を
   修正(移植版が本家より約2.1倍遅かった件の主因、修正済み)。**
   項目7でqssp180oldの`numerr=2`一致を確認した際、「なぜこの移植版は
   本家より約2.1倍遅いのか(実機1705秒 vs 移植版3557秒)」という
   追加の疑問が出たため、最初の5反復だけを`cProfile`で計測した
   (`pars["maxiter"]=5`で打ち切り、462秒)。判明した2つの問題点は
   いずれも計算結果に一切影響しない、純粋なPythonレベルのオーバー
   ヘッドだった:

   - **`_native.fwsolve()`/`bwsolve()`が`L_csc.indptr`/`.indices`/
     `.data`を毎回`np.ascontiguousarray(..., dtype=np.uintp)`で
     再変換していた。** `numeric_cholesky()`が返す`L_csc`
     (`scipy.sparse.csc_matrix`)は1回の外側反復のPCGループ全体で
     使い回される同一オブジェクトで、変わるのは右辺ベクトルだけ
     ―― にもかかわらず、`scipy.sparse.csc_matrix`はコンストラクタに
     `uintp`型のindicesを渡しても実際にはint32/int64に正規化して
     しまう(実際に確認済み)ため、`_as_index_array()`の`uintp`への
     再キャストが毎回のfwsolve/bwsolve呼び出しで律儀に発生していた。
     `qssp180old`は`nnz(L)~8.3e7`と非常に大きいため、
     `numpy.ascontiguousarray`だけでプロファイル対象462秒中130秒
     (28%)を占めていた。**修正**: 変換済み配列を`L_csc`オブジェクト
     自身の属性としてキャッシュする(`_cached_csc_solve_arrays()`)。
     `L_csc`は外側反復ごとに新しいオブジェクトが作られるので
     キャッシュが古くなることはなく、`deninfac.py`/`sdfactor.py`が
     `L_csc.data`/`.indices`を書き換えないことも確認済み。
   - **`_native.qblkmul()`(Lorentz錐ブロックごとのスカラー倍、
     `qblkmul.c`のmexFunctionに対応するCの独立関数が存在しないため
     NumPyへ直接移植していた箇所)が`for k in range(nblk)`という
     Pythonのforループでブロックを1つずつ処理していた。**
     `qssp180old`はLorentz錐ブロックが65,341個もあり、この関数の
     ループ本体だけでプロファイル対象462秒中49.7秒(334回呼び出し)を
     占めていた。**修正**: `np.repeat(mu, block_sizes) * d`で
     ベクトル化(各ブロックの`mu[k]`をブロック幅ぶん複製してから
     一括で掛け算する、同じ演算順序で同じ値を計算するだけの書き換え)。
     ベクトル化の過程で、一部の呼び出し元が「1ブロック分の長さぶん
     しか使わないのに、それより長い`d`を渡している」パターン
     (元のforループは黙って余分な末尾を無視していた)を見落として
     一度リグレッションを起こしたが、フルテストスイート(247件)が
     即座に検出したので`d = d[:span]`で明示的に揃えて解決。

   **効果**: `qssp180old`の最初の5反復のcProfileで**462秒→153秒
   (3.0倍高速化)**。`numpy.ascontiguousarray`の自己時間は130秒→
   4.2秒に、`qblkmul`はトップ項目から完全に消滅。残った最大のコストは
   `numeric_cholesky`(本物のCカーネルによるCholesky分解、75.5秒)
   ―― これが支配的になるのはむしろ健全(本家でも同様にここが
   支配的)。フルテストスイート(247件)、`pytest -m mini`
   (SDPLIB/DIMACS公開参照値との照合46問題)、`pytest -m extended`
   (同16問題、LP+SOCP比率の高い問題も含む)いずれも回帰なしで合格
   済み。qssp180oldの`numerr=2`一致自体(項目7)には影響なし
   (計算結果はbug-for-bug同一)。
9. **項目8の続き: 残っていた性能ボトルネックを調査 ―― ADA/DAt.qの
   疎行列再構築は本家自身の設計、`numeric_cholesky`のuintp/int64
   往復は一部のみ削減可能と判明。**

   修正後(項目8)の`qssp180old`最初5反復のcProfile(153秒)を
   さらに分析し、残っている主要コストの性質を切り分けた:

   - **ADA/DAt.qを毎反復scipy疎行列として一から組み立て直している
     コスト(`numpy.array`経由でscipy内部が8.8秒)は、本家SeDuMi自身の
     設計そのものであることを確認した。** `sum(K.s)==0`
     (`qssp180old`が通る経路)用の`vendor/sedumi-upstream/getada.m`
     (MEXではなく素の`.m`)を実際に読んだところ、`global ADA_sedumi_`
     を`ADA_sedumi_ = sparse([],[],[],m,m,nnz(ADA_sedumi_))`で**毎回
     空の疎行列として新規作成**し、`ADA_sedumi_ + DAt.q'*DAt.q`/
     `ADA_sedumi_ + Alq'*diag(sparse(scalingvector))*Alq`という疎×疎
     演算で組み立て直していた ―― まさにこの移植版の`getada.py`
     (`scipy.sparse`での同等の疎×疎演算)と同じパターン。
     `getDAtm.m`も`extractA`(MEX)で抽出した後
     `spdiags(d.q1,0,nq,nq) * DAt.q`(疎対角行列を毎回新規作成して
     掛ける)という、同じく都度組み立て直す実装だった。したがって
     このコストはポーティングの非効率ではなく**忠実な移植の結果**
     であり、修正の対象外と判断する。
     (余談: PSDブロックがある`sum(K.s)!=0`経路の`getada1.c`/
     `getada2.c`/`getada3.c`は対照的に、`getsymbada`が一度だけ確定
     させた疎パターンを使い回し、`ADA_sedumi_`という同じグローバル
     配列に**値だけ**書き込む「その場更新」設計になっている ――
     本家は「PSDありなら使い回す、PSDなしなら毎回作り直す」という
     非対称設計になっており、`qssp180old`が使う後者は本家の時点で
     最適化されていない、ということのようだ。)
   - **`numeric_cholesky`内の`Lir.astype(np.int64)`/
     `Ljc.astype(np.int64)`(項目8の`astype`18秒のほぼ全て)は、
     `scipy.sparse.csc_matrix`が索引配列をuintpでは保持せず必ず
     int32/int64に正規化してしまう(実測確認済み)ことに起因する、
     このコンテナ型を使う以上避けられないコスト**であり、項目8の
     見積もり(この`astype`往復を丸ごと削減できる)は誤りだった。
     実際に削減できたのは、`fwsolve()`/`bwsolve()`のキャッシュ
     (`_cached_csc_solve_arrays()`)が反復ごとの初回アクセス時に
     `L_csc.indptr`/`.indices`(int64)からuintpへ再変換していた分
     (1反復につき1回、計5回)だけだった。**修正**: `numeric_cholesky`
     が`L_csc`構築の直前にすでに持っているuintp版のLjc/Lirを、
     そのままキャッシュへ直接注入する(`L_csc._sedumipy_solve_cache`
     を`numeric_cholesky`自身が埋める)ことで、この最後の再変換を
     省略。**効果は153秒→138秒(約10%減)** ―― 当初見積もった16%には
     届かなかった。これ以上削るには`Lnum["L"]`をscipy疎行列から
     uintpをそのまま保持する自前の軽量構造体に置き換える必要がある
     (影響範囲がテストコード等にも及ぶ可能性がある、より大きな型
     変更のリファクタリングになるため、今回は見送り)。

   フルテストスイート(247件)、`pytest -m mini`(46問題)で回帰なし
   を確認済み。

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
