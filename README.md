# CESM2_PGWDEF

本實驗旨在分析區域海溫暖化與熱帶區域地表變遷對氣候的影響，採用地球系統模式 CESM2.3 進行實驗。

## 實驗設計

本實驗包括下列四組模擬：

- CTL：控制組。除了地表型態改為 1960 年的狀態外，所有條件均維持 B1850 compset 設定。
- DEF：地表變遷實驗。海洋大陸陸地區域地表型態改為 2020 年的狀態，而所有條件均維持 B1850 compset 設定。
- IO_DEF：地表變遷加海溫暖化實驗。海洋大陸陸地區域地表型態改為 2020 年的狀態，北印度洋海溫改為 pacemaker 驅動的暖化海溫條件，而所有條件均維持 B1850 compset 設定。
- IOWP_DEF：地表變遷加海溫暖化實驗。海洋大陸陸地區域地表型態改為 2020 年的狀態，北印度洋與西太平洋海溫改為 pacemaker 驅動的暖化海溫條件，而所有條件均維持 B1850 compset 設定。

## 附錄

### 地表型態數據

在本實驗中，CTL 與 DEF 使用 B1850 compset 模擬。為了使四組模擬均使用 1960 年的地表型態，本實驗使用了 BHIST compset 的地表數據，計算後建立了新的理想化地表數據。地表型態數據依照下列方式建立：
1. 讀取 flanduse 數據：landuse.timeseries_0.9x1.25_SSP5-8.5_78pfts_CMIP6_simyr1850-2100_c190214.nc，該數據為年資料。擷取其中的變數 PCT_NAT_PFT（以下簡稱 PFT）。
2. 擷取 1850 年、1960 年、2020 年的 PFT 資料。計算 1960-1850 年的差異值與 2020-1960 年的差異值。
3. 讀取 fsurdat 數據：B1850 需要使用 fsurdat 數據。因此，讀取 surfdata_0.9x1.25_hist_78pfts_CMIP6_simyr1850_c190214.nc，該數據為月資料（共 12 個月，季節循環）。擷取其中的變數 PCT_NAT_PFT。將該 PFT 加上前一步驟計算出來的 1960-1850 年差異值，作為 CTL 的地表型態數據。另外，將 PFT 加上 1960-1850 年差異值、並在海洋大陸陸地區域加上 2020-1960 年差異值，作為 DEF 的地表型態數據。

### 暖化數據

在本實驗中，區域海溫暖化使用了 Pacemaker 方式進行模擬。對於暖化實驗的模擬，一組「warming footprint」暖化海溫數據會被疊加至目標區域海域。warming footprint 數據依照下列方式建立：
1. 讀取 AMIP 實驗 TS 數據：AMIP 的海溫使用 ERSSTv5、海冰使用 HadISST1，該組實驗的 TS 數據可視為歷史觀測數據。在該數據中，擷取 1960 年至 2014 年的 TS 數據。
2. 每個網格點的趨勢值（單位為 K/year，使用回歸分析的最小平方法計算）。
