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
<td><img src="https://user-images.githubusercontent.com/109857341/194307518-df53ed6d-47fb-4a8c-a005-064abdea8af3.png" width="320px"></td>
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
<td><img src="https://user-images.githubusercontent.com/109857341/194307248-1eb3bf34-bc92-4e1e-8275-79690664ef3e.png" width="320px"></td>
</tr>
</table>

* 円柱周りの流れ（Re=100，層流）

https://user-images.githubusercontent.com/109857341/195817642-dc11b2af-c81f-4ab7-9215-0dbf3b18e614.mp4
<table>
<tr>
<td>Cx</td>
<td>Cy</td>
</tr>
<tr>
<td><img src="https://user-images.githubusercontent.com/109857341/195794865-2a9b0cf8-9825-48bc-81aa-85fd03c4bf30.png" width="320px"></td>
<td><img src="https://user-images.githubusercontent.com/109857341/195794938-69f6a37a-250e-42c5-a596-59c255b23292.png" width="320px"></td>
</tr>
</table>

* 動く円柱周りの流れ (層流)

https://user-images.githubusercontent.com/109857341/184158102-c2e7cd08-8ab5-4cbc-9b9f-9c97f7118309.mp4

* 後ろ向きステップの流れ（Re=5500，乱流）

https://user-images.githubusercontent.com/109857341/195793711-509993e8-f4a5-4960-bbc7-af8fb2dc5010.mp4

<img src="https://user-images.githubusercontent.com/109857341/195793783-89527b68-6543-470c-b5e6-5a5561d9a0c7.png" width="640px">


* 車周りの流れ（現実の物性値，乱流）

https://user-images.githubusercontent.com/109857341/194307950-79486366-6146-4bf7-a8bc-0f12c478207e.mp4

## ライセンス
[BSD-3-Clause license](LICENSE)
