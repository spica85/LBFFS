# LBFFS [[English](README.md)/Japanese]
LBFFSは格子ボルツマン法に基づく流体シミュレーションソフトです．C++とOpenCLで書かれており，単一のGPUで動作します．

## 特徴
* D3Q19離散速度モデル
* Recursive regularized衝突モデル
* Smagorinskyモデルに基づくLarge eddy simulation（壁近傍での乱流粘性の減衰を考慮）
* Half-way Bounce-Back法による格子点上での壁境界条件
* Filippova & Hanel’s Interpolated Bounce-Back法による格子点外での壁境界条件（曲面などへの対応）
* STLファイルによる壁境界の設定
* 反射波を抑制する流出境界条件 [Geier et al., Comput. Math. Appl. (2015), Appendix F]
* 流出境界からの反射波を抑制するスポンジゾーン
* 埋め込み境界法

## Google Colaboratoryでの実行スクリプト  
[Sample](runScriptOnColab.ipynb)

## テストケース
* ポアズイユ流れ (Re=100，層流)
<table>
<tr>
<td>Velocity distribution</td>
<td>u profile</td>
</tr>
<tr>
<td><img src="https://user-images.githubusercontent.com/109857341/180640617-7e83c0b4-61df-4ed4-ac4f-39554b86affe.png" width="320px"></td>
<td><img src="https://user-images.githubusercontent.com/109857341/180640633-b6779f8d-1921-493f-b64f-876f08a873d8.png" width="320px"></td>
</tr>
</table>

* 天井駆動のキャビティ流れ (Re=100，層流)
<table>
<tr>
<td>Velocity vector</td>
<td>u profile at x=0.5m</td>
</tr>
<tr>
<td><img src="https://user-images.githubusercontent.com/109857341/180638527-6905b752-ebff-4695-a5c2-aacec47b16ac.png" width="320px"></td>
<td><img src="https://user-images.githubusercontent.com/109857341/182006837-1e144cbc-5c16-4bc1-aefe-eba9ca35f386.png" width="320px"></td>
</tr>
</table>

* 円柱周りの流れ（Re=100，層流）

https://user-images.githubusercontent.com/109857341/180644337-b0e62fda-41a7-487d-9cee-98e37b96f939.mp4
<table>
<tr>
<td>Cx</td>
<td>Cy</td>
</tr>
<tr>
<td><img src="https://user-images.githubusercontent.com/109857341/180644297-db37a9b0-177e-4a2a-8390-9ada3c0c96dd.png" width="320px"></td>
<td><img src="https://user-images.githubusercontent.com/109857341/180644308-5bb19345-ec3a-4d10-8ad9-4452e091d878.png" width="320px"></td>
</tr>
</table>

* 動く円柱周りの流れ (層流)

https://user-images.githubusercontent.com/109857341/184156620-08c4a85d-8176-46f2-8cad-d64d2287c577.mp4

* 後ろ向きステップの流れ（Re=5500，乱流）

https://user-images.githubusercontent.com/109857341/180644458-212d29d3-9d87-4b73-b8d2-fd086c1d4b44.mp4

<img src="https://user-images.githubusercontent.com/109857341/180644496-94171507-2454-4ed6-b495-355ab656610b.png" width="640px">


* 車周りの流れ（現実の物性値，乱流）

https://user-images.githubusercontent.com/109857341/180644599-89a6945f-214d-449f-8b31-b7e3ab75fc98.mp4

## ライセンス
[BSD-3-Clause license](LICENSE)
