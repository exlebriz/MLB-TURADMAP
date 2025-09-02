# ============================================================================
# TÜRKİYE RADYONÜKLİD AKTİVİTE HARİTA & SORGU ARACI
# Version 5.1 — Düzeltilmiş ve Geliştirilmiş Versiyon
# Author: Cafer Mert YEŞİLKANAT 
# ============================================================================
Sys.setlocale()
Sys.setlocale(category = "LC_ALL", locale = "tr_TR.UTF-8")
### oncelikle calisma alanini tamimlayalim
# -------------------- GEREKLİ PAKETLER --------------------------------------
suppressPackageStartupMessages({
  library(shiny)
  library(shinyWidgets)
  library(bslib)
  library(DT)
  library(sf)
  library(terra)
  library(raster)
  library(leaflet)
  library(leaflet.extras)
  library(exactextractr)
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(gridExtra)
  library(plotly)
  library(corrplot)
  library(car)
  library(openxlsx)
})

options(shiny.maxRequestSize = 200*1024^2)   # 200 MB
LEAFLET_MAXBYTES <- 256 * 1024^2             # 256 MB
WEB_MERC <- "EPSG:3857"

# -------------------- DOSYA YOLLARI (kendine göre düzenle) ------------------
# -------------------- GÜNCEL DOSYA YOLLARI ------------------
PATH_RA_TIF <- "ra_3857.tif"
PATH_TH_TIF <- "th_3857.tif"
PATH_K_TIF  <- "k_3857.tif"
PATH_ID1_TIF <- "id1_3857.tif"
PATH_ID2_TIF <- "id2_3857.tif"
PATH_ILCE_MEANS_RDS <- "ilce_means.rds"

PATH_ILCE_SHP  <- "ilce_poligon.shp"
PATH_Ra_kmz <- "Ra226_XGBoost_RGB.kmz"
PATH_Th_kmz <- "Th232_XGBoost_RGB.kmz"
PATH_K_kmz  <- "K40_XGBoost_RGB.kmz"

# -------------------- RENK & EŞİK AYARLARI ----------------------------------
FORCED_MAX <- list(Ra=130, Th=108, K=1250)  
FORCED_MIN <- list(Ra=0, Th=0, K=0)         
PALETTE_FUN <- function(n) colorRampPalette(c("#0028ff","#25d9ff","#ffe14a","#ff8a1d","#ff2f2f"))(n)

# -------------------- YARDIMCI FONKSİYONLAR ---------------------------------
js_set_text <- HTML(
  "Shiny.addCustomMessageHandler('setText', function(x){
      var el = document.getElementById(x.id);
      if(el){ el.innerText = x.text; }
   });
   Shiny.addCustomMessageHandler('setButtonText', function(x){
      var btn = document.getElementById(x.id);
      if(btn && btn.tagName === 'A'){ 
        btn.innerText = x.text;
      }
   });"
)

safe_palette_numeric <- function(varname, domain_min, domain_max){
  cols <- PALETTE_FUN(256)
  leaflet::colorNumeric(palette = cols, domain = c(domain_min, domain_max), na.color = "transparent")
}

clip_to_max <- function(r, vmin, vmax){
  if (is.null(vmin) && is.null(vmax)) return(r)
  raster::calc(r, function(x){ 
    if (!is.null(vmin)) x[x < vmin] <- vmin
    if (!is.null(vmax)) x[x > vmax] <- vmax
    x 
  })
}

# Pearson Korelasyon Katsayısı hesaplama (CCC yerine)
calculate_pearson <- function(x, y) {
  if (length(x) != length(y)) return(NA)
  valid <- complete.cases(x, y)
  if (sum(valid) < 3) return(NA)
  
  x <- x[valid]
  y <- y[valid]
  
  pearson_r <- cor(x, y, method = "pearson")
  return(pearson_r)
}

# Mean Absolute Percentage Error (MAPE) hesaplama
# Mean Absolute Percentage Error (MAPE) hesaplama
# Mean Absolute Percentage Error (MAPE) hesaplama - OPTIMIZE EDİLMİŞ
calculate_mape <- function(actual, predicted) {
  if (length(actual) != length(predicted)) return(NA_real_)
  
  valid <- complete.cases(actual, predicted) & 
    is.finite(actual) & is.finite(predicted) & 
    actual != 0 & abs(actual) > 1e-10  # Sıfıra çok yakın değerleri filtrele
  
  if (sum(valid) < 2) return(NA_real_)
  
  actual <- actual[valid]
  predicted <- predicted[valid]
  
  # Vectorized hesaplama
  mape_vals <- abs((actual - predicted) / actual) * 100
  
  # Aşırı değerleri filtrele (outlier protection)
  mape_vals <- mape_vals[is.finite(mape_vals) & mape_vals < 1000]
  
  if (length(mape_vals) == 0) return(NA_real_)
  
  return(mean(mape_vals))
}

# Total Gamma Dose Rate hesaplama (nGy/h)
calculate_tgdr <- function(ra226, th232, k40) {
  # Standart katsayılar (nGy/h per Bq/kg)
  coeff_ra <- 0.462
  coeff_th <- 0.604 
  coeff_k <- 0.0417
  
  tgdr <- coeff_ra * ra226 + coeff_th * th232 + coeff_k * k40
  return(tgdr)
}

# İstatistiksel testler - Kolmogorov-Smirnov testi (Paired t-test yerine)
perform_statistical_tests <- function(observed, predicted, lang = "tr") {
  if (length(observed) < 3 || length(predicted) < 3) {
    return(data.frame(
      Test = c("Kolmogorov-Smirnov", "Wilcoxon signed-rank"),
      p_value = c(NA, NA),
      Result = c(
        if(lang=="tr") "Yetersiz veri" else "Insufficient data",
        if(lang=="tr") "Yetersiz veri" else "Insufficient data"
      ),
      stringsAsFactors = FALSE
    ))
  }
  
  # Kolmogorov-Smirnov testi (dağılım karşılaştırması)
  ks_test <- try(suppressWarnings(ks.test(observed, predicted)), silent = TRUE)
  ks_p_value <- ifelse(inherits(ks_test, "try-error"), NA, ks_test$p.value)
  ks_result <- ifelse(is.na(ks_p_value), 
                      if(lang=="tr") "Test başarısız" else "Test failed", 
                      ifelse(ks_p_value < 0.05, 
                             if(lang=="tr") "Dağılımlar farklı" else "Distributions different", 
                             if(lang=="tr") "Dağılımlar benzer" else "Distributions similar"))
  
  # Wilcoxon signed-rank test (eşleştirilmiş örnekler)
  wilcox_test <- try(wilcox.test(observed, predicted, paired = TRUE), silent = TRUE)
  w_p_value <- ifelse(inherits(wilcox_test, "try-error"), NA, wilcox_test$p.value)
  w_result <- ifelse(is.na(w_p_value), 
                     if(lang=="tr") "Test başarısız" else "Test failed", 
                     ifelse(w_p_value < 0.05, 
                            if(lang=="tr") "Anlamlı fark var" else "Significant difference", 
                            if(lang=="tr") "Anlamlı fark yok" else "No significant difference"))
  
  return(data.frame(
    Test = c("Kolmogorov-Smirnov", "Wilcoxon signed-rank"),
    p_value = round(c(ks_p_value, w_p_value), 4),
    Result = c(ks_result, w_result),
    stringsAsFactors = FALSE
  ))
}

# Elips hesaplama fonksiyonu (daha dar elips için)
calculate_ellipse <- function(x, y, confidence = 0.95) {  # 0.95'den 0.75'e düşürüldü
  if (length(x) < 3 || length(y) < 3) return(NULL)
  
  valid <- complete.cases(x, y)
  x <- x[valid]
  y <- y[valid]
  
  if (length(x) < 3) return(NULL)
  
  # Kovaryans matrisi
  data_matrix <- cbind(x, y)
  center <- colMeans(data_matrix)
  cov_matrix <- cov(data_matrix)
  
  # Chi-square değeri (daha dar elips)
  chi_sq <- qchisq(confidence, df = 2)
  
  # Eigen decomposition
  eigen_result <- eigen(cov_matrix)
  eigenvalues <- eigen_result$values
  eigenvectors <- eigen_result$vectors
  
  # Elips parametreleri
  a <- sqrt(chi_sq * eigenvalues[1])  # büyük eksen
  b <- sqrt(chi_sq * eigenvalues[2])  # küçük eksen
  
  # Açı
  angle <- atan2(eigenvectors[2,1], eigenvectors[1,1])
  
  # Elips noktaları
  theta <- seq(0, 2*pi, length.out = 100)
  ellipse_x <- center[1] + a * cos(theta) * cos(angle) - b * sin(theta) * sin(angle)
  ellipse_y <- center[2] + a * cos(theta) * sin(angle) + b * sin(theta) * cos(angle)
  
  return(list(
    x = ellipse_x,
    y = ellipse_y,
    center = center,
    a = a,
    b = b,
    angle = angle,
    confidence = confidence
  ))
}

# Noktanın elips içinde olup olmadığını kontrol eden fonksiyon
point_in_ellipse <- function(px, py, ellipse_info) {
  if (is.null(ellipse_info)) return(rep(FALSE, length(px)))
  
  # Merkeze göre normalize et
  x_norm <- px - ellipse_info$center[1]
  y_norm <- py - ellipse_info$center[2]
  
  # Rotasyon uygula
  cos_angle <- cos(-ellipse_info$angle)
  sin_angle <- sin(-ellipse_info$angle)
  
  x_rot <- x_norm * cos_angle - y_norm * sin_angle
  y_rot <- x_norm * sin_angle + y_norm * cos_angle
  
  # Elips formülü: (x/a)^2 + (y/b)^2 <= 1
  return((x_rot/ellipse_info$a)^2 + (y_rot/ellipse_info$b)^2 <= 1)
}

# TR/EN sözlük - GENİŞLETİLMİŞ
dict <- list(
  tr = list(
    app_title="Yüksek Çözünürlüklü Türkiye Radyonüklid Aktivite Haritası ve Sorgu Aracı",
    tab_query="Koordinat Sorgu",
    tab_explore="İl/İlçe Gezgini",
    tab_eval="İlçe Bazında Değerlendirme",
    map_title="Harita",
    map_hint="Haritaya tıklayarak koordinat ekleyebilirsiniz.",
    layer="Katman:",
    kmz_download="KMZ İndir",
    input_title="Koordinat Girişi",
    input_mode="Giriş Yöntemi:",
    manual="Manuel", paste="Kopyala-Yapıştır", file="Dosya Yükle",
    lon="Boylam:", lat="Enlem:",
    paste_placeholder="35.1234567,39.1234567",
    add="Ekle", query="Sorgula", clear="Temizle",
    results="Sonuçlar", download_results="Sonuçları İndir",
    region_title="Bölge Seçimi", 
    province="İl:", district="İlçe:", 
    multi_district="Çoklu İlçe Seçimi:",
    show_province="İli Göster",
    radionuclide="Radyonüklid:", show_map="Haritayı Göster",
    dist_avgs="Seçilen Bölge Ortalamaları",
    upload_title="İlçe Bazında Değerlendirme",
    eval_file="CSV Dosyası Seç:", 
    has_truth="Gerçek değerler mevcut (Ra226, Th232, K40)",
    evaluate="İlçe Bazında Değerlendir", 
    metrics="İlçe Karşılaştırma Metrikleri",
    plot_title="İlçe Bazında Çapraz Doğrulama",
    eval_results="İlçe Değerlendirme Sonuçları", 
    download_eval="Tüm Sonuçları İndir",
    download_plots="Grafikleri İndir (PDF)",
    ellipse_stats="İhtimaliyet Elipsi İstatistikleri",
    inside_ellipse="Elips İçi",
    outside_ellipse="Elips Dışı",
    district_comparison="İlçe Karşılaştırması",
    predicted_avg="Tahmin Ortalaması",
    observed_avg="Gözlem Ortalaması",
    point_count="Nokta Sayısı",
    ellipse_confidence="İhtimaliyet Elipsi (%95 güven)",
    district_ranking="İlçe Tahmin Performansı",
    ranking_desc="En iyi tahmin edilen ilçelerden en kötü tahmin edilen ilçelere doğru",
    coord_map_title="Yüklenen Koordinatlar Haritası",
    coord_map_hint="Tabloda satır seçerek haritada vurgulamaları değiştirebilirsiniz.",
    performance_score="Performans Skoru",
    rank="Sıra",
    district_metrics_title="İlçe Karşılaştırma Metrikleri",
    point_metrics_title="Nokta Bazında Metrikler", 
    tgdr_metrics_title="TGDR Metrikleri (nGy/h)",
    stat_tests_title="İstatistiksel Testler",
    point_plot_title="Nokta Bazında Çapraz Doğrulama",
    tgdr_plot_title="TGDR Çapraz Doğrulama",
    ellipse_inside_title="İhtimaliyet Elipsi İçi (%95)",
    point_ellipse_inside_title="İhtimaliyet Elipsi İçi (%95)",
    tgdr_ellipse_inside_title="İhtimaliyet Elipsi İçi (%95)",
    all_data="Tüm Veri Seti",
    all_points="Tüm Noktalar",
    all_tgdr="Tüm Veri",
    district_averages="İlçe Ortalamaları",
    point_based="Nokta Bazlı",
    tgdr_based="TGDR",
    kmz_missing="Seçili katman için KMZ dosyası bulunamadı.",
    csv_need_lonlat="CSV dosyası 'lon' ve 'lat' sütunlarını içermelidir!",
    no_data_district="Bu bölge için veri bulunamadı.",
    district_not_found="İlçe bulunamadı.",
    insufficient_data="Değerlendirme için yeterli veri yok.",
    evaluation_complete="İlçe bazında değerlendirme tamamlandı.",
    # İstatistiksel test açıklamaları
    stat_tests_explanation="İstatistiksel Testler: Kolmogorov-Smirnov testi gözlenen ve tahmin edilen değerlerin dağılımlarının benzerliğini test eder. Wilcoxon testi ise eşleştirilmiş örneklerin medyan değerlerini karşılaştırır.",
    # Yeni metrikler
    pearson_r="Pearson r"
  ),
  en = list(
    app_title="High-Resolution Türkiye Radionuclide Activity Map and Query Tool",
    tab_query="Coordinate Query",
    tab_explore="Province/District Explorer",
    tab_eval="District-Based Evaluation",
    map_title="Map",
    map_hint="Click on the map to add coordinates.",
    layer="Layer:",
    kmz_download="Download KMZ",
    input_title="Coordinate Input",
    input_mode="Input Mode:",
    manual="Manual", paste="Paste", file="Upload",
    lon="Longitude:", lat="Latitude:",
    paste_placeholder="35.1234567,39.1234567",
    add="Add", query="Query", clear="Clear",
    results="Results", download_results="Download Results",
    region_title="Region Selection", 
    province="Province:", district="District:",
    multi_district="Multiple District Selection:",
    show_province="Show Province",
    radionuclide="Radionuclide:", show_map="Show Map",
    dist_avgs="Selected Region Averages",
    upload_title="District-Based Evaluation",
    eval_file="Select CSV:", 
    has_truth="Ground truth present (Ra226, Th232, K40)",
    evaluate="District-Based Evaluate", 
    metrics="District Comparison Metrics",
    plot_title="District-Based Cross Validation",
    eval_results="District Evaluation Results", 
    download_eval="Download All Results",
    download_plots="Download Plots",
    ellipse_stats="Probability Ellipse Statistics",
    inside_ellipse="Inside Ellipse",
    outside_ellipse="Outside Ellipse",
    district_comparison="District Comparison",
    predicted_avg="Predicted Average",
    observed_avg="Observed Average",
    point_count="Point Count",
    ellipse_confidence="Probability Ellipse (95% confidence)",
    district_ranking="District Prediction Performance",
    ranking_desc="From best predicted districts to worst predicted districts",
    coord_map_title="Uploaded Coordinates Map",
    coord_map_hint="Select rows in the table to change highlighting on the map.",
    performance_score="Performance Score",
    rank="Rank",
    district_metrics_title="District Comparison Metrics",
    point_metrics_title="Point-Based Metrics", 
    tgdr_metrics_title="TGDR Metrics (nGy/h)",
    stat_tests_title="Statistical Tests",
    point_plot_title="Point-Based Cross Validation",
    tgdr_plot_title="TGDR Cross Validation",
    ellipse_inside_title="Inside Probability Ellipse (95%)",
    point_ellipse_inside_title="Inside Probability Ellipse (95%)",
    tgdr_ellipse_inside_title="Inside Probability Ellipse (95%)",
    all_data="All Dataset",
    all_points="All Points",
    all_tgdr="All Data",
    district_averages="District Averages",
    point_based="Point-Based",
    tgdr_based="TGDR",
    kmz_missing="No KMZ found for the selected layer.",
    csv_need_lonlat="CSV must contain 'lon' and 'lat' columns!",
    no_data_district="No data found for this region.",
    district_not_found="District not found.",
    insufficient_data="Insufficient data for evaluation.",
    evaluation_complete="District-based evaluation completed.",
    # İstatistiksel test açıklamaları
    stat_tests_explanation="Statistical Tests: Kolmogorov-Smirnov test compares the distributions of observed and predicted values. Wilcoxon test compares median values of paired samples.",
    # Yeni metrikler
    pearson_r="Pearson r"
  )
)

trn <- function(lang, key) {
  result <- dict[[lang]][[key]]
  if (is.null(result)) return(key)  # Anahtar bulunamazsa anahtarı döndür
  return(result)
}

# -------------------- VERİ YÜKLEME (Optimize Edilmiş) ---------------------
message("Ön-işlenmiş, optimize veriler yükleniyor...")

# terra paketinin raster'larını doğrudan ve hızlıca yükle
# Projeksiyon ve diğer ağır işlemler zaten yapıldı!
r_ra_web <- terra::rast(PATH_RA_TIF)
r_th_web <- terra::rast(PATH_TH_TIF)
r_k_web  <- terra::rast(PATH_K_TIF)

# Leaflet'in ihtiyaç duyduğu 'raster' paketi formatına anlık dönüştürme (bu hızlı bir işlem)
r_ra_r <- raster::raster(r_ra_web)
r_th_r <- raster::raster(r_th_web)
r_k_r  <- raster::raster(r_k_web)

# İl/İlçe verisini yükle
message("İl/İlçe sınırları yükleniyor...")
ilce_sf <- sf::st_read(PATH_ILCE_SHP, quiet = TRUE) |>
  sf::st_make_valid() |>
  sf::st_transform("EPSG:4326")

province_lookup <- ilce_sf |> sf::st_drop_geometry() |> dplyr::select(ID_1, NAME_1) |> dplyr::distinct() |> dplyr::arrange(NAME_1)
district_lookup <- ilce_sf |> sf::st_drop_geometry() |> dplyr::select(ID_1, NAME_1, ID_2, NAME_2) |> dplyr::distinct() |> dplyr::arrange(NAME_1, NAME_2)

# Önceden hesaplanmış ilçe ortalamalarını doğrudan yükle
message("Önceden hesaplanmış ilçe ortalamaları yükleniyor...")
ilce_means <- readRDS(PATH_ILCE_MEANS_RDS)

# Önceden hazırlanmış ID raster'larını yükle
message("ID rasterları yükleniyor...")
ID1_raster <- terra::rast(PATH_ID1_TIF)
ID2_raster <- terra::rast(PATH_ID2_TIF)

message("Tüm optimize veriler başarıyla yüklendi!")
# -------------------- VERİ YÜKLEME SONU --------------------------------

# -------------------- UI -----------------------------------------------------
ui <- fluidPage(
  theme = bslib::bs_theme(
    version = 5,
    base_font = bslib::font_google("Inter"),
    heading_font = bslib::font_google("Inter"),
    primary = "#6c5ce7", secondary = "#00cec9", success = "#00b894",
    info = "#0984e3", warning = "#fdcb6e", danger = "#d63031",
    bg = "#0f1117", fg = "#e6e6e6"
  ),
  tags$head(
    tags$style(HTML("
      .main-header{
        background: linear-gradient(135deg,#6c5ce7 0%,#00cec9 100%);
        color: white; padding: 20px; margin-bottom: 16px; border-radius: 14px;
      }
      .card{ background:#151923; border:1px solid #222836; border-radius:14px; padding:16px; margin-bottom:16px; }
      .leaflet-container{ border-radius:14px; box-shadow: 0 10px 24px rgba(0,0,0,.25); }
      .btn { border-radius: 10px; margin: 2px; }
      .metrics-card { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); }
      .stats-card { 
    background: linear-gradient(135deg, #2980b9 0%, #2c3e50 100%); /* Daha koyu ve okunur bir fon */
    color: white; /* Tüm kart içindeki yazı rengi beyaz olsun */
  }
  .stats-card h4, .stats-card p {
    color: white !important; /* Başlık ve paragraf renklerini zorla beyaz yap */
  }
      .results-card { background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%); }
      .stat-explanation { background-color: #2c3e50; color: #ecf0f1; padding: 10px; border-radius: 8px; margin: 10px 0; font-size: 0.9em; }
    ")),
    tags$script(js_set_text)
  ),
  
  # ÜST BAR
  div(class="main-header",
      fluidRow(
        column(8, h2(id="main-title","Türkiye Radyonüklid Aktivite Harita ve Sorgu Aracı")),
        column(4, align="right",
               radioButtons("lang", NULL, choices=c("Türkçe"="tr","English"="en"), selected="tr", inline=TRUE)
        )
      )
  ),
  
  
  
  
  tabsetPanel(id="tabs", type="tabs",
              # TAB 1: KOORDINAT SORGU
              tabPanel(title=span(id="tab1-title","Koordinat Sorgu"), value="query",
                       br(),
                       
                       # ÜST SATIR: SOLDa KONTROLLER (5), SAĞDA HARİTA (7)
                       fluidRow(
                         # SOL PANEL: Koordinat girişi + Katman & KMZ (hepsi solda)
                         column(4,
                                div(class="card",
                                    h4(id="input-title","Koordinat Girişi"),
                                    radioButtons("input_mode", label=span(id="input-mode-label","Giriş Yöntemi:"), 
                                                 choices=c("Manuel"="manual","Kopyala-Yapıştır"="paste","Dosya Yükle"="file"), 
                                                 selected="manual", inline=TRUE),
                                    hr(),
                                    uiOutput("input_ui"),
                                    br(),
                                    actionButton("run_query", label=span(id="query-button","Sorgula"), class="btn btn-primary"),
                                    actionButton("clear_points", label=span(id="clear-button","Temizle"), class="btn btn-warning")
                                ),
                                br(),
                                div(class="card",
                                    # Katman ve KMZ kontrolü de solda olsun (Tab 2 mantığına paralel)
                                    h5(span(id="layer-label","Katman:")),
                                    selectInput("overlay_layer", label=NULL, 
                                                choices=c("Ra-226"="ra","Th-232"="th","K-40"="k"), selected="ra"),
                                    div(align="right",
                                        downloadButton("download_kmz",label=span(id="kmz-button","KMZ İndir"), class="btn btn-info btn-sm")
                                    )
                                )
                         ),
                         
                         # SAĞ PANEL: Harita
                         column(8,
                                div(class="card",
                                    h4(id="map-title","Harita"),
                                    p(id="map-hint","Haritaya tıklayarak koordinat ekleyebilirsiniz."),
                                    leafletOutput("map_query", height="520px")
                                )
                         )
                       ),
                       
                       # ALT SATIR: SONUÇ TABLOSU TAM GENİŞLİKTE (12)
                       br(),
                       fluidRow(
                         column(12,
                                div(class="card",
                                    h4(id="results-title","Sonuçlar"),
                                    DT::dataTableOutput("results_table"),
                                    br(),
                                    conditionalPanel(condition="output.show_results_download",
                                                     downloadButton("download_results",label=span(id="download-results-button","Sonuçları İndir"), class="btn btn-success")
                                    )
                                )
                         )
                       )
              ),
              
              # TAB 2: GELİŞTİRİLMİŞ İL/İLÇE GEZGİNİ
              tabPanel(title=span(id="tab2-title","İl/İlçe Gezgini"), value="explore",
                       br(),
                       fluidRow(
                         column(4,
                                div(class="card",
                                    h4(id="region-title","Bölge Seçimi"),
                                    
                                    # İl seçimi
                                    selectInput("selected_province", label=span(id="province-label","İl:"), choices=NULL),
                                    
                                    # Sadece il gösterme butonu
                                    actionButton("show_province_only", label=span(id="show-province-button","İli Göster"), 
                                                 class="btn btn-success w-100", style="margin-bottom: 15px;"),
                                    
                                    hr(),
                                    
                                    # Çoklu ilçe seçimi
                                    h5(id="multi-district-title","Çoklu İlçe Seçimi:"),
                                    checkboxGroupInput("selected_districts","", choices=NULL),
                                    
                                    br(),
                                    selectInput("selected_nuclide", label=span(id="radionuclide-label","Radyonüklid:"), 
                                                choices=c("Ra-226"="ra","Th-232"="th","K-40"="k"), 
                                                selected="ra"),
                                    
                                    actionButton("show_districts", label=span(id="show-districts-button","İlçeleri Göster"), 
                                                 class="btn btn-primary w-100"),
                                    
                                    br(), br(),
                                    uiOutput("region_stats")
                                )
                         ),
                         column(8,
                                div(class="card",
                                    leafletOutput("map_explore", height="650px")
                                )
                         )
                       )
              ),
              
              # TAB 3: YENİ İLÇE BAZINDA DEĞERLENDİRME
              tabPanel(title=span(id="tab3-title","İlçe Bazında Değerlendirme"), value="evaluation",
                       br(),
                       fluidRow(
                         column(4,
                                # Dosya yükleme
                                div(class="card",
                                    h4(id="upload-title","İlçe Bazında Değerlendirme"),
                                    fileInput("eval_file", label=span(id="eval-file-label","CSV Dosyası Seç:"), accept = c(".csv", ".txt")),
                                    
                                    # Dosya önizleme
                                    conditionalPanel(condition="output.show_csv_preview",
                                                     tags$hr(),
                                                     h6("Dosya Önizleme:"),
                                                     verbatimTextOutput("csv_preview", placeholder = TRUE)
                                    ),
                                    
                                    checkboxInput("has_truth", label=span(id="has-truth-label","Gerçek değerler mevcut (Ra226, Th232, K40)"), value=FALSE),
                                    br(),
                                    actionButton("run_district_eval", label=span(id="evaluate-button","İlçe Bazında Değerlendir"), 
                                                 class="btn btn-primary w-100"),
                                    br(), br(),
                                    downloadButton("download_eval_results", label=span(id="download-eval-button","Tüm Sonuçları İndir"), 
                                                   class="btn btn-success w-100"),
                                    br(),
                                    downloadButton("download_eval_plots", label=span(id="download-plots-button","Grafikleri İndir"), 
                                                   class="btn btn-info w-100"),
                                    br(),
                                    downloadButton("download_eval_plots_svg", label=span(id="download-plots-svg-button","Grafikleri İndir (SVG ZIP)"),
                                                   class="btn btn-secondary w-100")  # <-- YENİ BUTON
                                ),
                                
                                # İstatistiksel Testler Açıklaması
                                conditionalPanel(condition="input.has_truth && output.show_stat_tests",
                                                 div(class="card stat-explanation",
                                                     h6(id="stat-test-info-title","İstatistiksel Testler Hakkında:"),
                                                     p(id="stat-explanation", 
                                                       "Kolmogorov-Smirnov testi: Gözlenen ve tahmin edilen değerlerin dağılımlarının benzerliğini test eder. Wilcoxon testi: Eşleştirilmiş örneklerin medyan değerlerini karşılaştırır.")
                                                 )
                                ),
                                
                                # İlçe Bazlı Metrikler
                                conditionalPanel(condition="input.has_truth && output.show_district_metrics",
                                                 div(class="card metrics-card",
                                                     h4(id="district-metrics-title","İlçe Karşılaştırma Metrikleri"),
                                                     h6(id="all-data-label","Tüm Veri Seti:"),
                                                     tableOutput("district_metrics_table"),
                                                     h6(id="ellipse-inside-title","İhtimaliyet Elipsi İçi (%95):"),
                                                     tableOutput("district_metrics_inside_table")
                                                 )
                                ),
                                
                                # Nokta Bazlı Metrikler
                                conditionalPanel(condition="input.has_truth && output.show_point_metrics",
                                                 div(class="card results-card",
                                                     h4(id="point-metrics-title","Nokta Bazında Metrikler"),
                                                     h6(id="all-points-label","Tüm Noktalar:"),
                                                     tableOutput("point_metrics_table"),
                                                     h6(id="point-ellipse-inside-title","İhtimaliyet Elipsi İçi (%95):"),
                                                     tableOutput("point_metrics_inside_table")
                                                 )
                                ),
                                
                                # TGDR Metrikleri
                                conditionalPanel(condition="input.has_truth && output.show_tgdr_metrics",
                                                 div(class="card stats-card",
                                                     h4(id="tgdr-metrics-title","TGDR Metrikleri (nGy/h)"),
                                                     h6(id="all-tgdr-label","Tüm Veri:"),
                                                     tableOutput("tgdr_metrics_table"),
                                                     h6(id="tgdr-ellipse-inside-title","İhtimaliyet Elipsi İçi (%95):"),
                                                     tableOutput("tgdr_metrics_inside_table")
                                                 )
                                ),
                                
                                # İstatistiksel Testler
                                # İstatistiksel Testler (Karşılaştırmalı)
                                conditionalPanel(condition="input.has_truth && output.show_stat_tests",
                                                 div(class="card",
                                                     h4(id="stat-tests-title","İstatistiksel Testler"),
                                                     
                                                     # ---- İLÇE BAZINDA ----
                                                     h5(id="district-averages-label","İlçe Ortalamaları:"),
                                                     h6(strong(style="color: #4facfe;", "Tüm Veri Seti:")),
                                                     tableOutput("district_stat_tests"),
                                                     h6(strong(style="color: #f093fb;", "Elips İçi Veri Seti (%95):")),
                                                     tableOutput("district_stat_tests_inside"),
                                                     hr(),
                                                     
                                                     # ---- NOKTA BAZINDA ----
                                                     h5(id="point-based-label","Nokta Bazlı:"),
                                                     h6(strong(style="color: #4facfe;", "Tüm Veri Seti:")),
                                                     tableOutput("point_stat_tests"),
                                                     h6(strong(style="color: #f093fb;", "Elips İçi Veri Seti (%95):")),
                                                     tableOutput("point_stat_tests_inside"),
                                                     hr(),
                                                     
                                                     # ---- TGDR BAZINDA ----
                                                     h5(id="tgdr-based-label","TGDR:"),
                                                     h6(strong(style="color: #4facfe;", "Tüm Veri Seti:")),
                                                     tableOutput("tgdr_stat_tests"),
                                                     h6(strong(style="color: #f093fb;", "Elips İçi Veri Seti (%95):")),
                                                     tableOutput("tgdr_stat_tests_inside")
                                                 )
                                ),
                                
                                # Elips istatistikleri
                                conditionalPanel(condition="input.has_truth && output.show_ellipse_stats",
                                                 div(class="card",
                                                     h4(id="ellipse-title","İhtimaliyet Elipsi İstatistikleri (%95)"),
                                                     h6(id="district-based-label","İlçe Bazında:"),
                                                     tableOutput("ellipse_stats_table"),
                                                     h6(id="point-based-ellipse-label","Nokta Bazında:"),
                                                     tableOutput("point_ellipse_stats_table"),
                                                     h6(id="tgdr-ellipse-label","TGDR:"),
                                                     tableOutput("tgdr_ellipse_stats_table")
                                                 )
                                )
                         ),
                         column(8,
                                # İlçe Bazlı Çapraz doğrulama grafikleri
                                conditionalPanel(condition="input.has_truth && output.show_district_plots",
                                                 div(class="card results-card",
                                                     h4(id="plot-title","İlçe Bazında Çapraz Doğrulama"),
                                                     plotOutput("district_validation_plots", height="500px")
                                                 )
                                ),
                                
                                # Nokta Bazlı Çapraz doğrulama grafikleri
                                conditionalPanel(condition="input.has_truth && output.show_point_plots",
                                                 div(class="card metrics-card",
                                                     h4(id="point-plot-title","Nokta Bazında Çapraz Doğrulama"),
                                                     plotOutput("point_validation_plots", height="500px")
                                                 )
                                ),
                                
                                # TGDR Çapraz doğrulama grafikleri
                                conditionalPanel(condition="input.has_truth && output.show_tgdr_plots",
                                                 div(class="card stats-card",
                                                     h4(id="tgdr-plot-title","TGDR Çapraz Doğrulama"),
                                                     plotOutput("tgdr_validation_plots", height="300px")
                                                 )
                                ),
                                
                                # İlçe Performans Sıralaması
                                conditionalPanel(condition="input.has_truth && output.show_district_ranking",
                                                 div(class="card stats-card",
                                                     h4(id="ranking-title","İlçe Tahmin Performansı Sıralaması"),
                                                     p(id="ranking-desc","En iyi tahmin edilen ilçelerden en kötü tahmin edilen ilçelere doğru sıralama"),
                                                     DT::dataTableOutput("district_ranking_table")
                                                 )
                                ),
                                
                                # Sonuç tablosu
                                div(class="card",
                                    h4(id="eval-results-title","İlçe Değerlendirme Sonuçları"),
                                    DT::dataTableOutput("district_eval_table")
                                )
                         ),
                         
                         # HARITA BÖLÜMÜ
                         conditionalPanel(condition="output.show_coord_map",
                                          column(12,
                                                 br(),
                                                 div(class="card",
                                                     h4(id="coord-map-title","Yüklenen Koordinatlar Haritası"),
                                                     p(id="coord-map-hint","Tabloda satır seçerek haritada vurgulamaları değiştirebilirsiniz."),
                                                     leafletOutput("coord_eval_map", height="600px")
                                                 )
                                          )
                         )
                       )
              )
  )
)

# -------------------- SERVER -------------------------------------------------
server <- function(input, output, session){
  
  # Dil durumu
  lang <- reactiveVal("tr")
  
  # Dil değişimi
  observeEvent(input$lang, {
    lang(input$lang)
  })
  
  # DİL GÜNCELLEMELERİ - KAPSAMLI VE DOĞRU
  observe({
    L <- lang()
    
    # Ana başlıklar ve tab isimleri
    session$sendCustomMessage("setText", list(id="main-title", text = trn(L, "app_title")))
    session$sendCustomMessage("setText", list(id="tab1-title", text = trn(L, "tab_query")))
    session$sendCustomMessage("setText", list(id="tab2-title", text = trn(L, "tab_explore")))
    session$sendCustomMessage("setText", list(id="tab3-title", text = trn(L, "tab_eval")))
    
    # Tab 1 elementleri
    session$sendCustomMessage("setText", list(id="map-title", text = trn(L, "map_title")))
    session$sendCustomMessage("setText", list(id="map-hint", text = trn(L, "map_hint")))
    session$sendCustomMessage("setText", list(id="input-title", text = trn(L, "input_title")))
    session$sendCustomMessage("setText", list(id="results-title", text = trn(L, "results")))
    session$sendCustomMessage("setText", list(id="layer-label", text = trn(L, "layer")))
    session$sendCustomMessage("setText", list(id="input-mode-label", text = trn(L, "input_mode")))
    
    # Tab 2 elementleri
    session$sendCustomMessage("setText", list(id="region-title", text = trn(L, "region_title")))
    session$sendCustomMessage("setText", list(id="multi-district-title", text = trn(L, "multi_district")))
    session$sendCustomMessage("setText", list(id="province-label", text = trn(L, "province")))
    session$sendCustomMessage("setText", list(id="radionuclide-label", text = trn(L, "radionuclide")))
    
    # Tab 3 elementleri
    session$sendCustomMessage("setText", list(id="upload-title", text = trn(L, "upload_title")))
    session$sendCustomMessage("setText", list(id="eval-file-label", text = trn(L, "eval_file")))
    session$sendCustomMessage("setText", list(id="has-truth-label", text = trn(L, "has_truth")))
    session$sendCustomMessage("setText", list(id="district-metrics-title", text = trn(L, "district_metrics_title")))
    session$sendCustomMessage("setText", list(id="point-metrics-title", text = trn(L, "point_metrics_title")))
    session$sendCustomMessage("setText", list(id="tgdr-metrics-title", text = trn(L, "tgdr_metrics_title")))
    session$sendCustomMessage("setText", list(id="stat-tests-title", text = trn(L, "stat_tests_title")))
    session$sendCustomMessage("setText", list(id="plot-title", text = trn(L, "plot_title")))
    session$sendCustomMessage("setText", list(id="point-plot-title", text = trn(L, "point_plot_title")))
    session$sendCustomMessage("setText", list(id="tgdr-plot-title", text = trn(L, "tgdr_plot_title")))
    session$sendCustomMessage("setText", list(id="eval-results-title", text = trn(L, "eval_results")))
    session$sendCustomMessage("setText", list(id="ellipse-title", text = trn(L, "ellipse_stats")))
    session$sendCustomMessage("setText", list(id="ranking-title", text = trn(L, "district_ranking")))
    session$sendCustomMessage("setText", list(id="ranking-desc", text = trn(L, "ranking_desc")))
    session$sendCustomMessage("setText", list(id="coord-map-title", text = trn(L, "coord_map_title")))
    session$sendCustomMessage("setText", list(id="coord-map-hint", text = trn(L, "coord_map_hint")))
    
    # Alt başlıklar
    session$sendCustomMessage("setText", list(id="ellipse-inside-title", text = trn(L, "ellipse_inside_title")))
    session$sendCustomMessage("setText", list(id="point-ellipse-inside-title", text = trn(L, "point_ellipse_inside_title")))
    session$sendCustomMessage("setText", list(id="tgdr-ellipse-inside-title", text = trn(L, "tgdr_ellipse_inside_title")))
    session$sendCustomMessage("setText", list(id="all-data-label", text = trn(L, "all_data")))
    session$sendCustomMessage("setText", list(id="all-points-label", text = trn(L, "all_points")))
    session$sendCustomMessage("setText", list(id="all-tgdr-label", text = trn(L, "all_tgdr")))
    session$sendCustomMessage("setText", list(id="district-averages-label", text = trn(L, "district_averages")))
    session$sendCustomMessage("setText", list(id="point-based-label", text = trn(L, "point_based")))
    session$sendCustomMessage("setText", list(id="tgdr-based-label", text = trn(L, "tgdr_based")))
    session$sendCustomMessage("setText", list(id="district-based-label", text = trn(L, "district_averages")))
    session$sendCustomMessage("setText", list(id="point-based-ellipse-label", text = trn(L, "point_based")))
    session$sendCustomMessage("setText", list(id="tgdr-ellipse-label", text = "TGDR"))
    
    # İstatistiksel test açıklaması
    session$sendCustomMessage("setText", list(id="stat-explanation", text = trn(L, "stat_tests_explanation")))
    session$sendCustomMessage("setText", list(id="stat-test-info-title", 
                                              text = if(L=="tr") "İstatistiksel Testler Hakkında:" else "About Statistical Tests:"))
    
    # Butonlar - setButtonText kullanarak
    session$sendCustomMessage("setButtonText", list(id="download_kmz", text = trn(L, "kmz_download")))
    session$sendCustomMessage("setButtonText", list(id="run_query", text = trn(L, "query")))
    session$sendCustomMessage("setButtonText", list(id="clear_points", text = trn(L, "clear")))
    session$sendCustomMessage("setButtonText", list(id="download_results", text = trn(L, "download_results")))
    session$sendCustomMessage("setButtonText", list(id="show_province_only", text = trn(L, "show_province")))
    session$sendCustomMessage("setButtonText", list(id="show_districts", text = trn(L, "show_map")))
    session$sendCustomMessage("setButtonText", list(id="run_district_eval", text = trn(L, "evaluate")))
    session$sendCustomMessage("setButtonText", list(id="download_eval_results", text = trn(L, "download_eval")))
    session$sendCustomMessage("setButtonText", list(id="download_eval_plots", text = trn(L, "download_plots")))
    
    # Input güncellemeleri
    updateRadioButtons(session, "input_mode", label = trn(L,"input_mode"),
                       choices = setNames(c("manual","paste","file"), c(trn(L,"manual"),trn(L,"paste"),trn(L,"file"))))
  })
  
  # Reaktif değerler
  rv <- reactiveValues(
    points = data.frame(lon = numeric(), lat = numeric()),
    next_id = 1,                    # <-- yeni sayaç
    query_results = NULL,
    district_eval_results = NULL,
    # İlçe bazlı metrikler
    district_metrics = NULL,
    district_metrics_inside = NULL,
    district_validation_plots = NULL,
    district_ellipse_stats = NULL,
    district_stat_tests = NULL,
    district_stat_tests_inside = NULL,
    # Nokta bazlı metrikler  
    point_metrics = NULL,
    point_metrics_inside = NULL,
    point_validation_plots = NULL,
    point_ellipse_stats = NULL,
    point_stat_tests = NULL,
    point_stat_tests_inside = NULL,
    # TGDR metrikler
    tgdr_metrics = NULL,
    tgdr_metrics_inside = NULL, 
    tgdr_validation_plots = NULL,
    tgdr_ellipse_stats = NULL,
    tgdr_stat_tests = NULL,
    tgdr_stat_tests_inside = NULL,
    # Diğer - SADECE BİR KEZ
    district_ranking = NULL,
    eval_coordinates = NULL,
    selected_rows = NULL
  )
  add_point_to_map <- function(lon, lat) {
    id_now <- rv$next_id
    rv$next_id <- rv$next_id + 1
    
    new_point <- data.frame(
      id  = id_now,
      lon = round(lon, 7),
      lat = round(lat, 7)
    )
    rv$points <- rbind(rv$points, new_point)
    
    # Kırmızı daire + kalıcı numara etiketi
    leafletProxy("map_query") |>
      addCircleMarkers(
        lng = new_point$lon, lat = new_point$lat,
        radius = 6, color = "#ff2f2f", fillOpacity = 0.85,
        layerId = paste0("pt_", id_now)
      ) |>
      addLabelOnlyMarkers(
        lng = new_point$lon, lat = new_point$lat,
        label = as.character(id_now),
        labelOptions = labelOptions(
          noHide = TRUE, direction = "top", textsize = "12px",
          style = list(
            "color" = "white", "font-weight" = "bold",
            "background" = "rgba(0,0,0,0.35)",
            "padding" = "2px 4px", "border-radius" = "6px"
          )
        ),
        layerId = paste0("lbl_", id_now)
      )
  }
  
  # ---- HARİTALAR ----
  output$map_query <- renderLeaflet({
    leaflet(options = leafletOptions(zoomSnap=0, zoomDelta=0.25)) |>
      addProviderTiles(providers$OpenStreetMap) |>
      setView(lng=35, lat=39, zoom=6) |>
      addPolygons(data=ilce_sf, color="gray", weight=1, fillOpacity=0, opacity=0.6, group="Districts")|>
      addScaleBar(position = "bottomleft",
                  options = list(metric = TRUE, imperial = FALSE, maxWidth = 200, updateWhenIdle = TRUE))
  })
  
  output$map_explore <- renderLeaflet({
    leaflet(options = leafletOptions(zoomSnap=0, zoomDelta=0.25)) |>
      addProviderTiles(providers$OpenStreetMap) |>
      setView(lng=35, lat=39, zoom=6) |>
      addPolygons(data=ilce_sf, color="gray", weight=1, fillOpacity=0, opacity=0.6, group="Districts")|>
      addScaleBar(position = "bottomleft",
                  options = list(metric = TRUE, imperial = FALSE, maxWidth = 200, updateWhenIdle = TRUE))
  })
  
  # Grid overlay ekleme fonksiyonu
  add_overlay <- function(mapId, layer_key){
    lyr <- switch(layer_key, ra=r_ra_r, th=r_th_r, k=r_k_r)
    ttl <- switch(layer_key, ra="Ra-226", th="Th-232", k="K-40")
    vmin <- switch(layer_key, ra=FORCED_MIN$Ra, th=FORCED_MIN$Th, k=FORCED_MIN$K)
    vmax <- switch(layer_key, ra=FORCED_MAX$Ra, th=FORCED_MAX$Th, k=FORCED_MAX$K)
    
    lyr <- clip_to_max(lyr, vmin, vmax)
    pal <- safe_palette_numeric(ttl, vmin, vmax)
    
    leafletProxy(mapId) |>
      clearImages() |>
      clearControls() |>
      addRasterImage(lyr, colors=pal, opacity=0.9, project=FALSE, maxBytes=LEAFLET_MAXBYTES) |>
      addLegend(pal=pal, values=c(vmin, vmax), title=paste(ttl,"(Bq/kg)"))|>
      addScaleBar(position = "bottomleft",
                  options = list(metric = TRUE, imperial = FALSE, maxWidth = 200, updateWhenIdle = TRUE))
  }
  
  # Başlangıç overlay
  observe({
    req(input$overlay_layer)
    add_overlay("map_query", input$overlay_layer)
  })
  
  observeEvent(input$overlay_layer, {
    add_overlay("map_query", input$overlay_layer)
  })
  
  # Harita tıklaması
  observeEvent(input$map_query_click, {
    click <- input$map_query_click
    add_point_to_map(click$lng, click$lat)
  })
  
  
  # ---- TAB 1 İŞLEVLERİ ----
  output$input_ui <- renderUI({
    req(input$input_mode)
    L <- lang()
    
    if (input$input_mode == "manual") {
      tagList(
        fluidRow(
          column(6, numericInput("man_lon", trn(L,"lon"), value=35, step=1e-7)),
          column(6, numericInput("man_lat", trn(L,"lat"), value=39, step=1e-7))
        ),
        actionButton("add_manual", trn(L,"add"), class="btn btn-info")
      )
    } else if (input$input_mode == "paste") {
      tagList(
        textAreaInput("paste_text", trn(L,"paste"),
                      rows=5, placeholder = trn(L,"paste_placeholder")),
        actionButton("add_paste", trn(L,"add"), class="btn btn-info")
      )
    } else {
      tagList(
        p(if (L=="tr") "CSV dosyası lon,lat sütunları içermelidir" else "CSV must include lon,lat columns"),
        fileInput("coords_file", if (L=="tr") "Koordinat Dosyası (.csv)" else "Coordinate File (.csv)",
                  accept=c(".csv",".txt")),
        actionButton("add_file", trn(L,"add"), class="btn btn-info")
      )
    }
  })
  
  observeEvent(input$add_manual, {
    req(input$man_lon, input$man_lat)
    add_point_to_map(input$man_lon, input$man_lat)
  })
  
  observeEvent(input$add_paste, {
    req(input$paste_text)
    lines <- strsplit(input$paste_text, "\n")[[1]]
    for (line in lines) {
      coords <- as.numeric(strsplit(trimws(line), ",")[[1]])
      if (length(coords)==2 && all(!is.na(coords))) {
        add_point_to_map(coords[1], coords[2])
      }
    }
  })
  
  
  observeEvent(input$add_file, {
    req(input$coords_file)
    df <- read.csv(input$coords_file$datapath)
    if (!all(c("lon","lat") %in% names(df))) {
      showNotification(trn(lang(),"csv_need_lonlat"), type="error"); return()
    }
    for (i in seq_len(nrow(df))) {
      add_point_to_map(df$lon[i], df$lat[i])
    }
  })
  
  
  observeEvent(input$clear_points, {
    rv$points <- data.frame(lon=numeric(), lat=numeric())
    rv$next_id <- 1                 # <-- sayaç reset
    rv$query_results <- NULL
    leafletProxy("map_query") |> clearMarkers()
  })
  
  # Sorgulama
  observeEvent(input$run_query, {
    req(nrow(rv$points)>0)
    withProgress(message = if (lang()=="tr") "Sorgulama yapılıyor..." else "Query running...", {
      pts_sf <- sf::st_as_sf(rv$points, coords=c("lon","lat"), crs="EPSG:4326")
      pts_3857 <- sf::st_transform(pts_sf, WEB_MERC)
      pts_vect <- terra::vect(pts_3857)
      
      ra_vals <- terra::extract(r_ra_web, pts_vect)[,2]
      th_vals <- terra::extract(r_th_web, pts_vect)[,2]
      k_vals  <- terra::extract(r_k_web,  pts_vect)[,2]
      
      id1 <- terra::extract(ID1_raster, pts_vect)[,2]
      id2 <- terra::extract(ID2_raster, pts_vect)[,2]
      
      results <- data.frame(
        ID = rv$points$id, 
        Longitude = rv$points$lon, Latitude = rv$points$lat,
        Ra226 = round(ra_vals,2), Th232 = round(th_vals,2), K40 = round(k_vals,2),
        Province = province_lookup$NAME_1[match(id1, province_lookup$ID_1)],
        District = district_lookup$NAME_2[match(id2, district_lookup$ID_2)]
      )
      
      results$Ra_Mean <- NA_real_; results$Th_Mean <- NA_real_; results$K_Mean <- NA_real_
      for (i in seq_len(nrow(results))) {
        if (!is.na(results$District[i])) {
          d <- ilce_means[ilce_means$NAME_2 == results$District[i], ]
          if (nrow(d)>0) {
            results$Ra_Mean[i] <- round(d$Ra_mean[1], 2)
            results$Th_Mean[i] <- round(d$Th_mean[1], 2)
            results$K_Mean[i]  <- round(d$K_mean[1], 2)
          }
        }
      }
      
      if (lang()=="tr") names(results) <- c("ID","Boylam","Enlem","Ra226","Th232","K40","İl","İlçe","Ra_Ort","Th_Ort","K_Ort")
      
      rv$query_results <- results
    })
  })
  
  output$results_table <- DT::renderDataTable({
    req(rv$query_results)
    DT::datatable(rv$query_results, options=list(pageLength=10))
  })
  
  output$show_results_download <- reactive({ 
    !is.null(rv$query_results) && nrow(rv$query_results)>0 
  })
  outputOptions(output, "show_results_download", suspendWhenHidden = FALSE)
  
  output$download_results <- downloadHandler(
    filename = function(){ paste0(if (lang()=="tr") "sonuclar_" else "results_", Sys.Date(), ".csv") },
    content = function(file){ write.csv(rv$query_results, file, row.names = FALSE) }
  )
  
  # KMZ indir
  output$download_kmz <- downloadHandler(
    filename = function() {
      nm <- switch(input$overlay_layer, ra="Ra226.kmz", th="Th232.kmz", k="K40.kmz")
      nm
    },
    content = function(file) {
      kmz_file <- switch(input$overlay_layer, ra=PATH_Ra_kmz, th=PATH_Th_kmz, k=PATH_K_kmz)
      if (file.exists(kmz_file)) {
        file.copy(kmz_file, file)
      } else { 
        showNotification(trn(lang(),"kmz_missing"), type="warning")
        writeLines("KMZ not available.", con=file) 
      }
    }
  )
  
  # ---- TAB 2: İL/İLÇE GEZGİNİ ----
  observe({
    provinces <- sort(unique(province_lookup$NAME_1))
    updateSelectInput(session, "selected_province", choices = provinces)
  })
  
  observeEvent(input$selected_province, {
    req(input$selected_province)
    districts <- district_lookup |> 
      dplyr::filter(NAME_1 == input$selected_province) |> 
      dplyr::pull(NAME_2) |> 
      sort()
    
    district_choices <- setNames(districts, districts)
    updateCheckboxGroupInput(session, "selected_districts", choices = district_choices)
  })
  
  # Sadece ili göster
  observeEvent(input$show_province_only, {
    req(input$selected_province, input$selected_nuclide)
    withProgress(message = if (lang()=="tr") "İl haritası yükleniyor..." else "Loading province map...", {
      
      # İl sınırını al
      province_poly <- ilce_sf |> 
        dplyr::filter(NAME_1 == input$selected_province)
      
      if (nrow(province_poly) == 0) {
        showNotification(trn(lang(),"district_not_found"), type="error")
        return()
      }
      
      # İl bazında birleştir
      province_union <- province_poly |> 
        sf::st_union() |> 
        sf::st_sf(geometry = _, NAME_1 = input$selected_province)
      
      # 3857'de kırp ve maskele
      province_union_3857 <- sf::st_transform(province_union, WEB_MERC)
      province_vect_3857 <- terra::vect(province_union_3857)
      
      raster_web <- switch(input$selected_nuclide, ra=r_ra_web, th=r_th_web, k=r_k_web)
      raster_crop <- terra::crop(raster_web, province_vect_3857)
      raster_mask <- terra::mask(raster_crop, province_vect_3857)
      raster_mask_r <- raster::raster(raster_mask)
      
      # Görselleştirme
      ttl  <- switch(input$selected_nuclide, ra="Ra-226", th="Th-232", k="K-40")
      vmin <- switch(input$selected_nuclide, ra=FORCED_MIN$Ra, th=FORCED_MIN$Th, k=FORCED_MIN$K)
      vmax <- switch(input$selected_nuclide, ra=FORCED_MAX$Ra, th=FORCED_MAX$Th, k=FORCED_MAX$K)
      raster_mask_r <- clip_to_max(raster_mask_r, vmin, vmax)
      pal <- safe_palette_numeric(ttl, vmin, vmax)
      
      # İl bazında istatistikler
      province_stats <- ilce_means |> 
        dplyr::filter(NAME_1 == input$selected_province) |>
        dplyr::summarise(
          Ra_mean = mean(Ra_mean, na.rm = TRUE),
          Th_mean = mean(Th_mean, na.rm = TRUE),
          K_mean  = mean(K_mean, na.rm = TRUE),
          district_count = n()
        )
      
      bbox <- sf::st_bbox(province_union)
      leafletProxy("map_explore") |>
        clearImages() |>
        clearShapes() |>
        clearControls() |>
        addRasterImage(raster_mask_r, colors=pal, opacity=0.9, project=FALSE, maxBytes=LEAFLET_MAXBYTES) |>
        addPolygons(data=province_union, color="white", weight=3, fillOpacity=0) |>
        addPolygons(data=province_poly, color="gray", weight=1, fillOpacity=0, opacity=0.5) |>
        addLegend(pal=pal, values=c(vmin, vmax), title=paste(ttl,"(Bq/kg)")) |>
        addScaleBar(position = "bottomleft",
                    options = list(metric = TRUE, imperial = FALSE, maxWidth = 200, updateWhenIdle = TRUE)) |>
        fitBounds(bbox[["xmin"]], bbox[["ymin"]], bbox[["xmax"]], bbox[["ymax"]])
      
      # İl istatistiklerini göster
      output$region_stats <- renderUI({
        tagList(
          h5(paste(input$selected_province, trn(lang(),"dist_avgs"))),
          HTML(sprintf("<p><strong>%s:</strong> %d</p>", 
                       if(lang()=="tr") "İlçe Sayısı" else "District Count", province_stats$district_count)),
          HTML(sprintf("<p><strong>Ra-226:</strong> %.2f Bq/kg</p>", province_stats$Ra_mean)),
          HTML(sprintf("<p><strong>Th-232:</strong> %.2f Bq/kg</p>", province_stats$Th_mean)),
          HTML(sprintf("<p><strong>K-40:</strong> %.2f Bq/kg</p>", province_stats$K_mean))
        )
      })
    })
  })
  
  # Seçili ilçeleri göster
  observeEvent(input$show_districts, {
    req(input$selected_province, input$selected_districts, input$selected_nuclide)
    req(length(input$selected_districts) > 0)
    
    withProgress(message = if (lang()=="tr") "İlçe haritaları yükleniyor..." else "Loading district maps...", {
      
      # Seçili ilçeleri al
      selected_district_polys <- ilce_sf |> 
        dplyr::filter(NAME_1 == input$selected_province, NAME_2 %in% input$selected_districts)
      
      if (nrow(selected_district_polys) == 0) {
        showNotification(trn(lang(),"district_not_found"), type="error")
        return()
      }
      
      # Birleştir
      districts_union <- selected_district_polys |> sf::st_union()
      
      # 3857'de işle
      districts_union_3857 <- sf::st_transform(districts_union, WEB_MERC) |> sf::st_sf(geometry = _)
      districts_vect_3857 <- terra::vect(districts_union_3857)
      
      raster_web <- switch(input$selected_nuclide, ra=r_ra_web, th=r_th_web, k=r_k_web)
      raster_crop <- terra::crop(raster_web, districts_vect_3857)
      raster_mask <- terra::mask(raster_crop, districts_vect_3857)
      raster_mask_r <- raster::raster(raster_mask)
      
      # Görselleştirme
      ttl  <- switch(input$selected_nuclide, ra="Ra-226", th="Th-232", k="K-40")
      vmin <- switch(input$selected_nuclide, ra=FORCED_MIN$Ra, th=FORCED_MIN$Th, k=FORCED_MIN$K)
      vmax <- switch(input$selected_nuclide, ra=FORCED_MAX$Ra, th=FORCED_MAX$Th, k=FORCED_MAX$K)
      raster_mask_r <- clip_to_max(raster_mask_r, vmin, vmax)
      pal <- safe_palette_numeric(ttl, vmin, vmax)
      
      # Seçili ilçeler için istatistikler
      districts_stats <- ilce_means |> 
        dplyr::filter(NAME_1 == input$selected_province, NAME_2 %in% input$selected_districts)
      
      bbox <- sf::st_bbox(selected_district_polys)
      leafletProxy("map_explore") |>
        clearImages() |>
        clearShapes() |>
        clearControls() |>
        addRasterImage(raster_mask_r, colors=pal, opacity=0.9, project=FALSE, maxBytes=LEAFLET_MAXBYTES) |>
        addPolygons(data=selected_district_polys, color="white", weight=2, fillOpacity=0) |>
        addPolygons(data=ilce_sf, color="gray", weight=1, fillOpacity=0, opacity=0.3) |>
        addLegend(pal=pal, values=c(vmin, vmax), title=paste(ttl,"(Bq/kg)")) |>
        addScaleBar(position = "bottomleft",
                    options = list(metric = TRUE, imperial = FALSE, maxWidth = 200, updateWhenIdle = TRUE)) |>
        fitBounds(bbox[["xmin"]], bbox[["ymin"]], bbox[["xmax"]], bbox[["ymax"]])
      
      # İlçe istatistiklerini göster
      output$region_stats <- renderUI({
        stats_list <- tagList(h5(trn(lang(),"dist_avgs")))
        
        for (i in seq_len(nrow(districts_stats))) {
          district_name <- districts_stats$NAME_2[i]
          ra_val <- districts_stats$Ra_mean[i]
          th_val <- districts_stats$Th_mean[i]
          k_val  <- districts_stats$K_mean[i]
          
          stats_list <- tagList(stats_list,
                                h6(strong(district_name)),
                                HTML(sprintf("<p>Ra-226: %.2f | Th-232: %.2f | K-40: %.2f Bq/kg</p>", 
                                             ra_val, th_val, k_val)),
                                hr()
          )
        }
        
        stats_list
      })
    })
  })
  
  # ---- CSV ÖNIZLEME ----
  output$csv_preview <- renderText({
    req(input$eval_file)
    
    tryCatch({
      # Dosyayı ham olarak oku
      raw_lines <- readLines(input$eval_file$datapath, n = 5, warn = FALSE)
      
      # CSV olarak oku
      df_preview <- read.csv(input$eval_file$datapath, nrows = 3, stringsAsFactors = FALSE)
      
      # Koordinat tiplerini kontrol et
      lon_samples <- head(df_preview[["lon"]], 3) 
      lat_samples <- head(df_preview[["lat"]], 3)
      
      preview_text <- paste0(
        "=== ", if(lang()=="tr") "DOSYA BİLGİLERİ" else "FILE INFORMATION", " ===\n",
        if(lang()=="tr") "Ham satırlar:" else "Raw lines:", "\n", 
        paste(raw_lines[1:min(3, length(raw_lines))], collapse = "\n"), "\n\n",
        if(lang()=="tr") "Sütunlar:" else "Columns:", " ", paste(names(df_preview), collapse = ", "), "\n",
        if(lang()=="tr") "Koordinat örnekleri:" else "Coordinate samples:", "\n",
        "  lon: ", paste(lon_samples, collapse = ", "), "\n",
        "  lat: ", paste(lat_samples, collapse = ", "), "\n",
        if(lang()=="tr") "Koordinat sınıfları:" else "Coordinate classes:", "\n",
        "  lon class: ", class(df_preview[["lon"]]), "\n",
        "  lat class: ", class(df_preview[["lat"]]), "\n"
      )
      
      # Gerçek değer kontrolü
      if (input$has_truth) {
        truth_cols <- c("Ra226", "Th232", "K40")
        available_truth <- truth_cols[truth_cols %in% names(df_preview)]
        preview_text <- paste0(preview_text, 
                               if(lang()=="tr") "Gerçek değer sütunları:" else "Ground truth columns:", " ", 
                               paste(available_truth, collapse = ", "), "\n")
      }
      
      return(preview_text)
    }, error = function(e) {
      return(paste(if(lang()=="tr") "Dosya okuma hatası:" else "File reading error:", e$message))
    })
  })
  
  output$show_csv_preview <- reactive({
    !is.null(input$eval_file)
  })
  outputOptions(output, "show_csv_preview", suspendWhenHidden = FALSE)
  
  # ---- TAB 3: İLÇE BAZINDA DEĞERLENDİRME ----
  observeEvent(input$run_district_eval, {
    req(input$eval_file)
    
    withProgress(message = if (lang()=="tr") "İlçe bazında değerlendirme yapılıyor..." else "District-based evaluation running...", 
                 value = 0, {
                   
                   # CSV oku
                   incProgress(0.1, detail = if (lang()=="tr") "Dosya okunuyor..." else "Reading file...")
                   
                   tryCatch({
                     df <- read.csv(input$eval_file$datapath, stringsAsFactors = FALSE)
                   }, error = function(e) {
                     showNotification(paste("CSV okuma hatası:", e$message), type="error")
                     return(NULL)
                   })
                   
                   if (is.null(df)) return()
                   
                   # Sütun isimlerini kontrol et
                   if (!all(c("lon","lat") %in% names(df))) { 
                     available_cols <- paste(names(df), collapse = ", ")
                     showNotification(paste(trn(lang(),"csv_need_lonlat"), if(lang()=="tr") "Mevcut sütunlar:" else "Available columns:", available_cols), type="error")
                     return() 
                   }
                   
                   # Koordinat sütunlarını temizle ve kontrol et
                   incProgress(0.15, detail = if (lang()=="tr") "Koordinatlar kontrol ediliyor..." else "Checking coordinates...")
                   
                   # Koordinatları numeric'e dönüştür ve NA kontrolü yap
                   df$lon_original <- df$lon
                   df$lat_original <- df$lat
                   
                   # Önce string temizliği yap
                   df$lon <- gsub(",", ".", as.character(df$lon))  # virgülü noktaya çevir
                   df$lat <- gsub(",", ".", as.character(df$lat))  # virgülü noktaya çevir
                   
                   # Boşlukları temizle
                   df$lon <- trimws(df$lon)
                   df$lat <- trimws(df$lat)
                   
                   # Numeric'e dönüştür
                   df$lon <- suppressWarnings(as.numeric(df$lon))
                   df$lat <- suppressWarnings(as.numeric(df$lat))
                   
                   # Orijinal satır sayısı
                   original_rows <- nrow(df)
                   
                   # NA koordinatları temizle
                   df <- df[complete.cases(df[c("lon", "lat")]), ]
                   clean_rows <- nrow(df)
                   
                   if (clean_rows == 0) {
                     showNotification(if (lang()=="tr") "Geçerli koordinat bulunamadı! Koordinatlar sayısal formatta olmalı." 
                                      else "No valid coordinates found! Coordinates must be in numeric format.", type="error")
                     return()
                   }
                   
                   # Koordinat aralığını kontrol et (Türkiye sınırları)
                   valid_coords <- df$lon >= 25 & df$lon <= 45 & df$lat >= 35 & df$lat <= 43
                   df <- df[valid_coords, ]
                   
                   if (nrow(df) == 0) {
                     showNotification(if (lang()=="tr") "Koordinatlar Türkiye sınırları dışında! (Boylam: 25-45, Enlem: 35-43)" 
                                      else "Coordinates outside Turkey boundaries! (Lon: 25-45, Lat: 35-43)", type="error")
                     return()
                   }
                   
                   # Kullanıcıya bilgi ver
                   if (clean_rows < original_rows) {
                     showNotification(paste(if (lang()=="tr") "Temizlenen veriler:" else "Data cleaning:", 
                                            original_rows - clean_rows, 
                                            if (lang()=="tr") "geçersiz koordinat silindi" else "invalid coordinates removed"), 
                                      type="warning", duration = 5)
                   }
                   
                   # Her koordinatın hangi ilçede olduğunu bul
                   incProgress(0.2, detail = if (lang()=="tr") "İlçeler belirleniyor..." else "Identifying districts...")
                   
                   tryCatch({
                     pts_sf <- sf::st_as_sf(df[,c("lon","lat")], coords=c("lon","lat"), crs="EPSG:4326")
                   }, error = function(e) {
                     showNotification(paste("Koordinat dönüşüm hatası:", e$message), type="error")
                     return(NULL)
                   })
                   
                   if (is.null(pts_sf)) return()
                   
                   # Spatial join ile ilçe bilgisini al
                   incProgress(0.25, detail = if (lang()=="tr") "İlçe eşleştirmeleri yapılıyor..." else "Matching districts...")
                   
                   tryCatch({
                     pts_with_districts <- sf::st_join(pts_sf, ilce_sf)
                   }, error = function(e) {
                     showNotification(paste("İlçe eşleştirme hatası:", e$message), type="error")
                     return(NULL)
                   })
                   
                   if (is.null(pts_with_districts)) return()
                   
                   # Orijinal data'ya ilçe bilgilerini ekle
                   df$ID_1 <- pts_with_districts$ID_1
                   df$NAME_1 <- pts_with_districts$NAME_1
                   df$ID_2 <- pts_with_districts$ID_2
                   df$NAME_2 <- pts_with_districts$NAME_2
                   
                   # Grid'den tahmin değerlerini al
                   incProgress(0.3, detail = if (lang()=="tr") "Tahmin değerleri alınıyor..." else "Extracting predictions...")
                   
                   tryCatch({
                     pts_3857 <- sf::st_transform(pts_sf, WEB_MERC)
                     pts_vect <- terra::vect(pts_3857)
                     
                     df$Ra226_grid <- terra::extract(r_ra_web, pts_vect)[,2]
                     df$Th232_grid <- terra::extract(r_th_web, pts_vect)[,2]
                     df$K40_grid   <- terra::extract(r_k_web,  pts_vect)[,2]
                   }, error = function(e) {
                     showNotification(paste("Grid tahmin değerleri alınırken hata:", e$message), type="error")
                     return(NULL)
                   })
                   
                   # İlçe sınırları içinde olan noktaları filtrele
                   df <- df[complete.cases(df[c("NAME_1", "NAME_2")]), ]
                   
                   # Grid değerleri içinde NA olanları da temizle
                   df <- df[complete.cases(df[c("Ra226_grid", "Th232_grid", "K40_grid")]), ]
                   
                   if (nrow(df) == 0) {
                     showNotification(if (lang()=="tr") "İlçe sınırları içinde hiç koordinat bulunamadı!" 
                                      else "No coordinates found within district boundaries!", type="error")
                     return()
                   }
                   
                   # Bilgi mesajı
                   showNotification(paste(if (lang()=="tr") "İşlenen koordinat sayısı:" else "Processed coordinates:", nrow(df)), 
                                    type="default", duration = 3)
                   
                   incProgress(0.5, detail = if (lang()=="tr") "İlçe ortalamaları hesaplanıyor..." else "Calculating district averages...")
                   
                   # İlçe bazında gruplama
                   if (input$has_truth && all(c("Ra226","Th232","K40") %in% names(df))) {
                     
                     # Gerçek değerleri de kontrol et ve temizle
                     if ("Ra226" %in% names(df)) {
                       df$Ra226 <- gsub(",", ".", as.character(df$Ra226))
                       df$Ra226 <- trimws(df$Ra226)
                       df$Ra226 <- suppressWarnings(as.numeric(df$Ra226))
                     }
                     
                     if ("Th232" %in% names(df)) {
                       df$Th232 <- gsub(",", ".", as.character(df$Th232))
                       df$Th232 <- trimws(df$Th232)
                       df$Th232 <- suppressWarnings(as.numeric(df$Th232))
                     }
                     
                     if ("K40" %in% names(df)) {
                       df$K40 <- gsub(",", ".", as.character(df$K40))
                       df$K40 <- trimws(df$K40)
                       df$K40 <- suppressWarnings(as.numeric(df$K40))
                     }
                     
                     # Gerçek değerlerde NA olanları temizle
                     df <- df[complete.cases(df[c("Ra226", "Th232", "K40")]), ]
                     
                     if (nrow(df) == 0) {
                       showNotification(if (lang()=="tr") "Gerçek değerlerde geçerli veri bulunamadı!" 
                                        else "No valid data found in ground truth values!", type="error")
                       return()
                     }
                     
                     # Her ilçe için gözlem ve tahmin ortalamalarını hesapla
                     district_summary <- df |>
                       group_by(ID_1, NAME_1, ID_2, NAME_2) |>
                       summarise(
                         # Gözlem ortalamaları
                         Ra226_observed_mean = round(mean(Ra226, na.rm = TRUE),2),
                         Th232_observed_mean = round(mean(Th232, na.rm = TRUE),2),
                         K40_observed_mean   = round(mean(K40, na.rm = TRUE),2),
                         
                         # Grid tahmin ortalamaları
                         Ra226_pixel_mean = round(mean(Ra226_grid, na.rm = TRUE),2),
                         Th232_pixel_mean = round(mean(Th232_grid, na.rm = TRUE),2),
                         K40_pixel_mean   = round(mean(K40_grid, na.rm = TRUE),2),
                         
                         # Nokta sayısı
                         point_count = n(),
                         .groups = 'drop'
                       )
                     
                     # İlçe tahmin ortalamalarını (genel harita ortalamaları) ekle
                     district_summary <- district_summary |>
                       left_join(ilce_means |> select(ID_2, Ra_mean, Th_mean, K_mean), by = "ID_2") |>
                       rename(
                         Ra226_map_mean = Ra_mean,
                         Th232_map_mean = Th_mean,
                         K40_map_mean   = K_mean
                       )
                     
                     # En az 2 nokta olan ilçeleri al
                     district_summary <- district_summary |> filter(point_count >= 1)
                     
                     if (nrow(district_summary) < 3) {
                       showNotification(trn(lang(),"insufficient_data"), type="error")
                       return()
                     }
                     
                     # === 1. İLÇE BAZINDA DEĞERLENDİRME ===
                     incProgress(0.7, detail = if (lang()=="tr") "İlçe bazlı metrikler hesaplanıyor..." else "Calculating district-based metrics...")
                     
                     # İlçe bazında performans metrikleri hesapla (Pearson r ile)
                     district_metrics <- data.frame(
                       Radionuclide = c("Ra-226", "Th-232", "K-40"),
                       RMSE = c(
                         sqrt(mean((district_summary$Ra226_observed_mean - district_summary$Ra226_map_mean)^2, na.rm = TRUE)),
                         sqrt(mean((district_summary$Th232_observed_mean - district_summary$Th232_map_mean)^2, na.rm = TRUE)),
                         sqrt(mean((district_summary$K40_observed_mean - district_summary$K40_map_mean)^2, na.rm = TRUE))
                       ),
                       MAPE = c(
                         calculate_mape(district_summary$Ra226_observed_mean, district_summary$Ra226_map_mean),
                         calculate_mape(district_summary$Th232_observed_mean, district_summary$Th232_map_mean),
                         calculate_mape(district_summary$K40_observed_mean, district_summary$K40_map_mean)
                       ),
                       MAE = c(
                         mean(abs(district_summary$Ra226_observed_mean - district_summary$Ra226_map_mean), na.rm = TRUE),
                         mean(abs(district_summary$Th232_observed_mean - district_summary$Th232_map_mean), na.rm = TRUE),
                         mean(abs(district_summary$K40_observed_mean - district_summary$K40_map_mean), na.rm = TRUE)
                       ),
                       Pearson_r = c(
                         calculate_pearson(district_summary$Ra226_observed_mean, district_summary$Ra226_map_mean),
                         calculate_pearson(district_summary$Th232_observed_mean, district_summary$Th232_map_mean),
                         calculate_pearson(district_summary$K40_observed_mean, district_summary$K40_map_mean)
                       )
                     )
                     
                     rv$district_metrics <- district_metrics
                     
                     # İlçe bazında elips hesaplamaları ve içi/dışı karşılaştırma
                     district_ellipse_ra <- calculate_ellipse(district_summary$Ra226_map_mean, district_summary$Ra226_observed_mean, 0.95)
                     district_ellipse_th <- calculate_ellipse(district_summary$Th232_map_mean, district_summary$Th232_observed_mean, 0.95)
                     district_ellipse_k <- calculate_ellipse(district_summary$K40_map_mean, district_summary$K40_observed_mean, 0.95)
                     
                     # Elips içi noktalar
                     district_inside_ra <- if(!is.null(district_ellipse_ra)) point_in_ellipse(district_summary$Ra226_map_mean, district_summary$Ra226_observed_mean, district_ellipse_ra) else rep(FALSE, nrow(district_summary))
                     district_inside_th <- if(!is.null(district_ellipse_th)) point_in_ellipse(district_summary$Th232_map_mean, district_summary$Th232_observed_mean, district_ellipse_th) else rep(FALSE, nrow(district_summary))
                     district_inside_k <- if(!is.null(district_ellipse_k)) point_in_ellipse(district_summary$K40_map_mean, district_summary$K40_observed_mean, district_ellipse_k) else rep(FALSE, nrow(district_summary))
                     
                     # İlçe bazında elips içi metrikler
                     district_metrics_inside <- data.frame(
                       Radionuclide = c("Ra-226", "Th-232", "K-40"),
                       RMSE = c(
                         if(sum(district_inside_ra) >= 2) sqrt(mean((district_summary$Ra226_observed_mean[district_inside_ra] - district_summary$Ra226_map_mean[district_inside_ra])^2, na.rm = TRUE)) else NA,
                         if(sum(district_inside_th) >= 2) sqrt(mean((district_summary$Th232_observed_mean[district_inside_th] - district_summary$Th232_map_mean[district_inside_th])^2, na.rm = TRUE)) else NA,
                         if(sum(district_inside_k) >= 2) sqrt(mean((district_summary$K40_observed_mean[district_inside_k] - district_summary$K40_map_mean[district_inside_k])^2, na.rm = TRUE)) else NA
                       ),
                       MAPE = c(
                         if(sum(district_inside_ra) >= 2) calculate_mape(district_summary$Ra226_observed_mean[district_inside_ra], district_summary$Ra226_map_mean[district_inside_ra]) else NA,
                         if(sum(district_inside_th) >= 2) calculate_mape(district_summary$Th232_observed_mean[district_inside_th], district_summary$Th232_map_mean[district_inside_th]) else NA,
                         if(sum(district_inside_k) >= 2) calculate_mape(district_summary$K40_observed_mean[district_inside_k], district_summary$K40_map_mean[district_inside_k]) else NA
                       ),
                       MAE = c(
                         if(sum(district_inside_ra) >= 2) mean(abs(district_summary$Ra226_observed_mean[district_inside_ra] - district_summary$Ra226_map_mean[district_inside_ra]), na.rm = TRUE) else NA,
                         if(sum(district_inside_th) >= 2) mean(abs(district_summary$Th232_observed_mean[district_inside_th] - district_summary$Th232_map_mean[district_inside_th]), na.rm = TRUE) else NA,
                         if(sum(district_inside_k) >= 2) mean(abs(district_summary$K40_observed_mean[district_inside_k] - district_summary$K40_map_mean[district_inside_k]), na.rm = TRUE) else NA
                       ),
                       Pearson_r = c(
                         if(sum(district_inside_ra) >= 2) calculate_pearson(district_summary$Ra226_observed_mean[district_inside_ra], district_summary$Ra226_map_mean[district_inside_ra]) else NA,
                         if(sum(district_inside_th) >= 2) calculate_pearson(district_summary$Th232_observed_mean[district_inside_th], district_summary$Th232_map_mean[district_inside_th]) else NA,
                         if(sum(district_inside_k) >= 2) calculate_pearson(district_summary$K40_observed_mean[district_inside_k], district_summary$K40_map_mean[district_inside_k]) else NA
                       )
                     )
                     
                     rv$district_metrics_inside <- district_metrics_inside
                     
                     # İlçe bazında elips istatistikleri
                     rv$district_ellipse_stats <- data.frame(
                       Radionuclide = c("Ra-226", "Th-232", "K-40"),
                       Inside_Ellipse = c(sum(district_inside_ra), sum(district_inside_th), sum(district_inside_k)),
                       Total_Points = c(length(district_inside_ra), length(district_inside_th), length(district_inside_k)),
                       Percentage_Inside = round(c(
                         sum(district_inside_ra) / length(district_inside_ra) * 100,
                         sum(district_inside_th) / length(district_inside_th) * 100,
                         sum(district_inside_k) / length(district_inside_k) * 100
                       ), 1)
                     )
                     
                     # İlçe bazında istatistiksel testler (Kolmogorov-Smirnov ile)
                     # İlçe bazında istatistiksel testler (Kolmogorov-Smirnov ile)
                     current_lang <- lang()
                     
                     # 1. Tüm Veri Seti ile Testler
                     rv$district_stat_tests <- rbind(
                       data.frame(Radionuclide = "Ra-226", perform_statistical_tests(district_summary$Ra226_observed_mean, district_summary$Ra226_map_mean, current_lang)),
                       data.frame(Radionuclide = "Th-232", perform_statistical_tests(district_summary$Th232_observed_mean, district_summary$Th232_map_mean, current_lang)),
                       data.frame(Radionuclide = "K-40", perform_statistical_tests(district_summary$K40_observed_mean, district_summary$K40_map_mean, current_lang))
                     )
                     
                     # 2. SADECE ELİPS İÇİ Veri Seti ile Testler
                     rv$district_stat_tests_inside <- rbind(
                       data.frame(Radionuclide = "Ra-226", perform_statistical_tests(
                         district_summary$Ra226_observed_mean[district_inside_ra], 
                         district_summary$Ra226_map_mean[district_inside_ra], 
                         current_lang)),
                       data.frame(Radionuclide = "Th-232", perform_statistical_tests(
                         district_summary$Th232_observed_mean[district_inside_th], 
                         district_summary$Th232_map_mean[district_inside_th], 
                         current_lang)),
                       data.frame(Radionuclide = "K-40", perform_statistical_tests(
                         district_summary$K40_observed_mean[district_inside_k], 
                         district_summary$K40_map_mean[district_inside_k], 
                         current_lang))
                     )
                     
                     # === 2. NOKTA BAZINDA DEĞERLENDİRME ===
                     incProgress(0.75, detail = if (lang()=="tr") "Nokta bazlı metrikler hesaplanıyor..." else "Calculating point-based metrics...")
                     
                     # Nokta bazında performans metrikleri (Pearson r ile)
                     point_metrics <- data.frame(
                       Radionuclide = c("Ra-226", "Th-232", "K-40"),
                       RMSE = c(
                         sqrt(mean((df$Ra226 - df$Ra226_grid)^2, na.rm = TRUE)),
                         sqrt(mean((df$Th232 - df$Th232_grid)^2, na.rm = TRUE)),
                         sqrt(mean((df$K40 - df$K40_grid)^2, na.rm = TRUE))
                       ),
                       MAPE = c(
                         calculate_mape(df$Ra226, df$Ra226_grid),
                         calculate_mape(df$Th232, df$Th232_grid),
                         calculate_mape(df$K40, df$K40_grid)
                       ),
                       MAE = c(
                         mean(abs(df$Ra226 - df$Ra226_grid), na.rm = TRUE),
                         mean(abs(df$Th232 - df$Th232_grid), na.rm = TRUE),
                         mean(abs(df$K40 - df$K40_grid), na.rm = TRUE)
                       ),
                       Pearson_r = c(
                         calculate_pearson(df$Ra226, df$Ra226_grid),
                         calculate_pearson(df$Th232, df$Th232_grid),
                         calculate_pearson(df$K40, df$K40_grid)
                       )
                     )
                     
                     rv$point_metrics <- point_metrics
                     
                     # Nokta bazında elips hesaplamaları
                     point_ellipse_ra <- calculate_ellipse(df$Ra226_grid, df$Ra226, 0.95)
                     point_ellipse_th <- calculate_ellipse(df$Th232_grid, df$Th232, 0.95)
                     point_ellipse_k <- calculate_ellipse(df$K40_grid, df$K40, 0.95)
                     
                     # Elips içi noktalar
                     point_inside_ra <- if(!is.null(point_ellipse_ra)) point_in_ellipse(df$Ra226_grid, df$Ra226, point_ellipse_ra) else rep(FALSE, nrow(df))
                     point_inside_th <- if(!is.null(point_ellipse_th)) point_in_ellipse(df$Th232_grid, df$Th232, point_ellipse_th) else rep(FALSE, nrow(df))
                     point_inside_k <- if(!is.null(point_ellipse_k)) point_in_ellipse(df$K40_grid, df$K40, point_ellipse_k) else rep(FALSE, nrow(df))
                     
                     # Nokta bazında elips içi metrikler
                     point_metrics_inside <- data.frame(
                       Radionuclide = c("Ra-226", "Th-232", "K-40"),
                       RMSE = c(
                         if(sum(point_inside_ra) >= 5) sqrt(mean((df$Ra226[point_inside_ra] - df$Ra226_grid[point_inside_ra])^2, na.rm = TRUE)) else NA,
                         if(sum(point_inside_th) >= 5) sqrt(mean((df$Th232[point_inside_th] - df$Th232_grid[point_inside_th])^2, na.rm = TRUE)) else NA,
                         if(sum(point_inside_k) >= 5) sqrt(mean((df$K40[point_inside_k] - df$K40_grid[point_inside_k])^2, na.rm = TRUE)) else NA
                       ),
                       MAPE = c(
                         if(sum(point_inside_ra) >= 5) calculate_mape(df$Ra226[point_inside_ra], df$Ra226_grid[point_inside_ra]) else NA,
                         if(sum(point_inside_th) >= 5) calculate_mape(df$Th232[point_inside_th], df$Th232_grid[point_inside_th]) else NA,
                         if(sum(point_inside_k) >= 5) calculate_mape(df$K40[point_inside_k], df$K40_grid[point_inside_k]) else NA
                       ),
                       MAE = c(
                         if(sum(point_inside_ra) >= 5) mean(abs(df$Ra226[point_inside_ra] - df$Ra226_grid[point_inside_ra]), na.rm = TRUE) else NA,
                         if(sum(point_inside_th) >= 5) mean(abs(df$Th232[point_inside_th] - df$Th232_grid[point_inside_th]), na.rm = TRUE) else NA,
                         if(sum(point_inside_k) >= 5) mean(abs(df$K40[point_inside_k] - df$K40_grid[point_inside_k]), na.rm = TRUE) else NA
                       ),
                       Pearson_r = c(
                         if(sum(point_inside_ra) >= 5) calculate_pearson(df$Ra226[point_inside_ra], df$Ra226_grid[point_inside_ra]) else NA,
                         if(sum(point_inside_th) >= 5) calculate_pearson(df$Th232[point_inside_th], df$Th232_grid[point_inside_th]) else NA,
                         if(sum(point_inside_k) >= 5) calculate_pearson(df$K40[point_inside_k], df$K40_grid[point_inside_k]) else NA
                       )
                     )
                     
                     rv$point_metrics_inside <- point_metrics_inside
                     
                     # Nokta bazında elips istatistikleri
                     rv$point_ellipse_stats <- data.frame(
                       Radionuclide = c("Ra-226", "Th-232", "K-40"),
                       Inside_Ellipse = c(sum(point_inside_ra), sum(point_inside_th), sum(point_inside_k)),
                       Total_Points = c(length(point_inside_ra), length(point_inside_th), length(point_inside_k)),
                       Percentage_Inside = round(c(
                         sum(point_inside_ra) / length(point_inside_ra) * 100,
                         sum(point_inside_th) / length(point_inside_th) * 100,
                         sum(point_inside_k) / length(point_inside_k) * 100
                       ), 1)
                     )
                     
                     # Nokta bazında istatistiksel testler (Kolmogorov-Smirnov ile)
                     # Nokta bazında istatistiksel testler (Kolmogorov-Smirnov ile)
                     
                     # 1. Tüm Veri Seti ile Testler
                     rv$point_stat_tests <- rbind(
                       data.frame(Radionuclide = "Ra-226", perform_statistical_tests(df$Ra226, df$Ra226_grid, current_lang)),
                       data.frame(Radionuclide = "Th-232", perform_statistical_tests(df$Th232, df$Th232_grid, current_lang)),
                       data.frame(Radionuclide = "K-40", perform_statistical_tests(df$K40, df$K40_grid, current_lang))
                     )
                     
                     # 2. SADECE ELİPS İÇİ Veri Seti ile Testler
                     rv$point_stat_tests_inside <- rbind(
                       data.frame(Radionuclide = "Ra-226", perform_statistical_tests(
                         df$Ra226[point_inside_ra], 
                         df$Ra226_grid[point_inside_ra], 
                         current_lang)),
                       data.frame(Radionuclide = "Th-232", perform_statistical_tests(
                         df$Th232[point_inside_th], 
                         df$Th232_grid[point_inside_th], 
                         current_lang)),
                       data.frame(Radionuclide = "K-40", perform_statistical_tests(
                         df$K40[point_inside_k], 
                         df$K40_grid[point_inside_k], 
                         current_lang))
                     )
                     
                     # === 3. TGDR DEĞERLENDİRMESİ ===
                     incProgress(0.85, detail = if (lang()=="tr") "TGDR metrikleri hesaplanıyor..." else "Calculating TGDR metrics...")
                     
                     # TGDR hesapla
                     df$TGDR_observed <- calculate_tgdr(df$Ra226, df$Th232, df$K40)
                     df$TGDR_predicted <- calculate_tgdr(df$Ra226_grid, df$Th232_grid, df$K40_grid)
                     
                     # TGDR performans metrikleri (Pearson r ile)
                     tgdr_metrics <- data.frame(
                       Metric = c("RMSE", "MAPE", "MAE", "Pearson r"),
                       Value = c(
                         sqrt(mean((df$TGDR_observed - df$TGDR_predicted)^2, na.rm = TRUE)),
                         calculate_mape(df$TGDR_observed, df$TGDR_predicted),
                         mean(abs(df$TGDR_observed - df$TGDR_predicted), na.rm = TRUE),
                         calculate_pearson(df$TGDR_observed, df$TGDR_predicted)
                       )
                     )
                     
                     rv$tgdr_metrics <- tgdr_metrics
                     
                     # TGDR elips hesaplaması
                     tgdr_ellipse <- calculate_ellipse(df$TGDR_predicted, df$TGDR_observed, 0.95)
                     
                     # Elips içi noktalar
                     tgdr_inside <- if(!is.null(tgdr_ellipse)) point_in_ellipse(df$TGDR_predicted, df$TGDR_observed, tgdr_ellipse) else rep(FALSE, nrow(df))
                     
                     # TGDR elips içi metrikler
                     tgdr_metrics_inside <- if(sum(tgdr_inside) >= 5) {
                       data.frame(
                         Metric = c("RMSE", "MAPE", "MAE", "Pearson r"),
                         Value = c(
                           sqrt(mean((df$TGDR_observed[tgdr_inside] - df$TGDR_predicted[tgdr_inside])^2, na.rm = TRUE)),
                           calculate_mape(df$TGDR_observed[tgdr_inside], df$TGDR_predicted[tgdr_inside]),
                           mean(abs(df$TGDR_observed[tgdr_inside] - df$TGDR_predicted[tgdr_inside]), na.rm = TRUE),
                           calculate_pearson(df$TGDR_observed[tgdr_inside], df$TGDR_predicted[tgdr_inside])
                         )
                       )
                     } else {
                       data.frame(
                         Metric = c("RMSE", "MAPE", "MAE", "Pearson r"),
                         Value = c(NA, NA, NA, NA)
                       )
                     }
                     
                     rv$tgdr_metrics_inside <- tgdr_metrics_inside
                     
                     # TGDR elips istatistikleri
                     rv$tgdr_ellipse_stats <- data.frame(
                       Metric = "TGDR",
                       Inside_Ellipse = sum(tgdr_inside),
                       Total_Points = length(tgdr_inside),
                       Percentage_Inside = round(sum(tgdr_inside) / length(tgdr_inside) * 100, 1)
                     )
                     
                     # TGDR istatistiksel testler (Kolmogorov-Smirnov ile)
                     # TGDR istatistiksel testler (Kolmogorov-Smirnov ile)
                     
                     # 1. Tüm Veri Seti ile Testler
                     rv$tgdr_stat_tests <- data.frame(
                       Metric = "TGDR",
                       perform_statistical_tests(df$TGDR_observed, df$TGDR_predicted, current_lang)
                     )
                     
                     # 2. SADECE ELİPS İÇİ Veri Seti ile Testler
                     rv$tgdr_stat_tests_inside <- data.frame(
                       Metric = "TGDR",
                       perform_statistical_tests(
                         df$TGDR_observed[tgdr_inside], 
                         df$TGDR_predicted[tgdr_inside], 
                         current_lang)
                     )
                     
                     # === 4. GRAFİKLER OLUŞTURMA ===
                     incProgress(0.9, detail = if (lang()=="tr") "Grafikler hazırlanıyor..." else "Preparing plots...")
                     
                     # İlçe bazında grafik fonksiyonu (kare görünüm)
                     create_district_plots <- function() {
                       par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3.5, 1), pty = "s")  # pty="s" kare yapar
                       
                       # Eksenleri geniş tutmak için range hesapla
                       ra_range <- range(c(district_summary$Ra226_map_mean, district_summary$Ra226_observed_mean), na.rm = TRUE)
                       ra_expand <- diff(ra_range) * 0.2
                       ra_xlim <- ra_ylim <- c(ra_range[1] - ra_expand, ra_range[2] + ra_expand)
                       
                       th_range <- range(c(district_summary$Th232_map_mean, district_summary$Th232_observed_mean), na.rm = TRUE)  
                       th_expand <- diff(th_range) * 0.2
                       th_xlim <- th_ylim <- c(th_range[1] - th_expand, th_range[2] + th_expand)
                       
                       k_range <- range(c(district_summary$K40_map_mean, district_summary$K40_observed_mean), na.rm = TRUE)
                       k_expand <- diff(k_range) * 0.2
                       k_xlim <- k_ylim <- c(k_range[1] - k_expand, k_range[2] + k_expand)
                       
                       # Ra-226 (kare)
                       plot(district_summary$Ra226_map_mean, district_summary$Ra226_observed_mean,
                            xlim = ra_xlim, ylim = ra_ylim, asp = 1,
                            xlab = if(lang()=="tr") "Harita Ortalaması (Bq/kg)" else "Map Average (Bq/kg)",
                            ylab = if(lang()=="tr") "Gözlem Ortalaması (Bq/kg)" else "Observed Average (Bq/kg)",
                            main = "Ra-226", pch = 19, cex = 1.2,
                            col = ifelse(district_inside_ra, "#e74c3c", "#3498db"))
                       
                       if (!is.null(district_ellipse_ra)) {
                         lines(district_ellipse_ra$x, district_ellipse_ra$y, col = "#f39c12", lwd = 2, lty = 2)
                       }
                       abline(0, 1, col = "#2c3e50", lwd = 2)
                       grid(col = "#bdc3c7", lty = 3)
                       
                       # Th-232 (kare)
                       plot(district_summary$Th232_map_mean, district_summary$Th232_observed_mean,
                            xlim = th_xlim, ylim = th_ylim, asp = 1,
                            xlab = if(lang()=="tr") "Harita Ortalaması (Bq/kg)" else "Map Average (Bq/kg)",
                            ylab = if(lang()=="tr") "Gözlem Ortalaması (Bq/kg)" else "Observed Average (Bq/kg)",
                            main = "Th-232", pch = 19, cex = 1.2,
                            col = ifelse(district_inside_th, "#e74c3c", "#27ae60"))
                       
                       if (!is.null(district_ellipse_th)) {
                         lines(district_ellipse_th$x, district_ellipse_th$y, col = "#f39c12", lwd = 2, lty = 2)
                       }
                       abline(0, 1, col = "#2c3e50", lwd = 2)
                       grid(col = "#bdc3c7", lty = 3)
                       
                       # K-40 (kare)
                       plot(district_summary$K40_map_mean, district_summary$K40_observed_mean,
                            xlim = k_xlim, ylim = k_ylim, asp = 1,
                            xlab = if(lang()=="tr") "Harita Ortalaması (Bq/kg)" else "Map Average (Bq/kg)",
                            ylab = if(lang()=="tr") "Gözlem Ortalaması (Bq/kg)" else "Observed Average (Bq/kg)",
                            main = "K-40", pch = 19, cex = 1.2,
                            col = ifelse(district_inside_k, "#e74c3c", "#9b59b6"))
                       
                       if (!is.null(district_ellipse_k)) {
                         lines(district_ellipse_k$x, district_ellipse_k$y, col = "#f39c12", lwd = 2, lty = 2)
                       }
                       abline(0, 1, col = "#2c3e50", lwd = 2)
                       grid(col = "#bdc3c7", lty = 3)
                     }
                     
                     # Nokta bazında grafik fonksiyonu (kare görünüm) 
                     create_point_plots <- function() {
                       par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3.5, 1), pty = "s")  # pty="s" kare yapar
                       
                       # Geniş aralıklar
                       ra_range <- range(c(df$Ra226_grid, df$Ra226), na.rm = TRUE)
                       ra_expand <- diff(ra_range) * 0.3
                       ra_xlim <- ra_ylim <- c(ra_range[1] - ra_expand, ra_range[2] + ra_expand)
                       
                       th_range <- range(c(df$Th232_grid, df$Th232), na.rm = TRUE)  
                       th_expand <- diff(th_range) * 0.3
                       th_xlim <- th_ylim <- c(th_range[1] - th_expand, th_range[2] + th_expand)
                       
                       k_range <- range(c(df$K40_grid, df$K40), na.rm = TRUE)
                       k_expand <- diff(k_range) * 0.3
                       k_xlim <- k_ylim <- c(k_range[1] - k_expand, k_range[2] + k_expand)
                       
                       # Ra-226 (kare)
                       plot(df$Ra226_grid, df$Ra226, xlim = ra_xlim, ylim = ra_ylim, asp = 1,
                            xlab = if(lang()=="tr") "Tahmin (Bq/kg)" else "Predicted (Bq/kg)",
                            ylab = if(lang()=="tr") "Gözlem (Bq/kg)" else "Observed (Bq/kg)",
                            main = "Ra-226", pch = 16, cex = 0.8, 
                            col = ifelse(point_inside_ra, "#e74c3c", "#3498db"))
                       
                       if (!is.null(point_ellipse_ra)) {
                         lines(point_ellipse_ra$x, point_ellipse_ra$y, col = "#f39c12", lwd = 2, lty = 2)
                       }
                       abline(0, 1, col = "#2c3e50", lwd = 2)
                       grid(col = "#bdc3c7", lty = 3)
                       
                       # Th-232 (kare)
                       plot(df$Th232_grid, df$Th232, xlim = th_xlim, ylim = th_ylim, asp = 1,
                            xlab = if(lang()=="tr") "Tahmin (Bq/kg)" else "Predicted (Bq/kg)",
                            ylab = if(lang()=="tr") "Gözlem (Bq/kg)" else "Observed (Bq/kg)",
                            main = "Th-232", pch = 16, cex = 0.8,
                            col = ifelse(point_inside_th, "#e74c3c", "#27ae60"))
                       
                       if (!is.null(point_ellipse_th)) {
                         lines(point_ellipse_th$x, point_ellipse_th$y, col = "#f39c12", lwd = 2, lty = 2)
                       }
                       abline(0, 1, col = "#2c3e50", lwd = 2)
                       grid(col = "#bdc3c7", lty = 3)
                       
                       # K-40 (kare)
                       plot(df$K40_grid, df$K40, xlim = k_xlim, ylim = k_ylim, asp = 1,
                            xlab = if(lang()=="tr") "Tahmin (Bq/kg)" else "Predicted (Bq/kg)",
                            ylab = if(lang()=="tr") "Gözlem (Bq/kg)" else "Observed (Bq/kg)",
                            main = "K-40", pch = 16, cex = 0.8,
                            col = ifelse(point_inside_k, "#e74c3c", "#9b59b6"))
                       
                       if (!is.null(point_ellipse_k)) {
                         lines(point_ellipse_k$x, point_ellipse_k$y, col = "#f39c12", lwd = 2, lty = 2)
                       }
                       abline(0, 1, col = "#2c3e50", lwd = 2)
                       grid(col = "#bdc3c7", lty = 3)
                     }
                     
                     # TGDR grafik fonksiyonu (kare görünüm)
                     create_tgdr_plot <- function() {
                       par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3.5, 2), pty = "s")  # pty="s" kare yapar
                       
                       # Geniş aralık
                       tgdr_range <- range(c(df$TGDR_predicted, df$TGDR_observed), na.rm = TRUE)
                       tgdr_expand <- diff(tgdr_range) * 0.3
                       tgdr_xlim <- tgdr_ylim <- c(tgdr_range[1] - tgdr_expand, tgdr_range[2] + tgdr_expand)
                       
                       plot(df$TGDR_predicted, df$TGDR_observed, xlim = tgdr_xlim, ylim = tgdr_ylim, asp = 1,
                            xlab = if(lang()=="tr") "Tahmin TGDR (nGy/h)" else "Predicted TGDR (nGy/h)",
                            ylab = if(lang()=="tr") "Gözlem TGDR (nGy/h)" else "Observed TGDR (nGy/h)",
                            main = if(lang()=="tr") "TGDR Çapraz Doğrulama" else "TGDR Cross Validation", 
                            pch = 16, cex = 1.0,
                            col = ifelse(tgdr_inside, "#e74c3c", "#8e44ad"))
                       
                       if (!is.null(tgdr_ellipse)) {
                         lines(tgdr_ellipse$x, tgdr_ellipse$y, col = "#f39c12", lwd = 3, lty = 2)
                       }
                       abline(0, 1, col = "#2c3e50", lwd = 2)
                       grid(col = "#bdc3c7", lty = 3)
                     }
                     
                     rv$district_validation_plots <- create_district_plots
                     rv$point_validation_plots <- create_point_plots  
                     rv$tgdr_validation_plots <- create_tgdr_plot
                     
                     # İlçe bazında performans sıralaması hesapla (MAPE bazında)
                     incProgress(0.95, detail = if (lang()=="tr") "Performans sıralaması hazırlanıyor..." else "Preparing performance ranking...")
                     
               
                     # İlçe bazında performans sıralaması hesapla (sadece MAPE bazında)
                     # İlçe bazında performans sıralaması hesapla (sadece MAPE bazında)
                     # İlçe bazında performans sıralaması hesapla (sadece MAPE bazında)
                     # Vectorized hesaplama - çok daha hızlı
                     district_performance <- district_summary %>%
                       mutate(
                         Ra226_MAPE = purrr::map2_dbl(Ra226_observed_mean, Ra226_map_mean, calculate_mape),
                         Th232_MAPE = purrr::map2_dbl(Th232_observed_mean, Th232_map_mean, calculate_mape),
                         K40_MAPE   = purrr::map2_dbl(K40_observed_mean, K40_map_mean, calculate_mape)
                       ) %>%
                       rowwise() %>%
                       mutate(Mean_MAPE = mean(c(Ra226_MAPE, Th232_MAPE, K40_MAPE), na.rm = TRUE)) %>%
                       ungroup() %>%
                       filter(!is.na(Mean_MAPE) & is.finite(Mean_MAPE)) %>%
                       arrange(Mean_MAPE) %>%
                       mutate(
                         Rank = row_number(),
                         Performance_Category = case_when(
                           Mean_MAPE <= quantile(Mean_MAPE, 0.33, na.rm = TRUE) ~ if(input$lang=="tr") "İyi" else "Good",
                           Mean_MAPE <= quantile(Mean_MAPE, 0.66, na.rm = TRUE) ~ if(input$lang=="tr") "Orta" else "Average", 
                           TRUE ~ if(input$lang=="tr") "Zayıf" else "Poor"
                         )
                       ) %>%
                       select(Rank, NAME_1, NAME_2, point_count, Mean_MAPE, Performance_Category, 
                              Ra226_MAPE, Th232_MAPE, K40_MAPE)
                     
                     # Güvenli atama
                     if(nrow(district_performance) > 0) {
                       rv$district_ranking <- district_performance
                     } else {
                       rv$district_ranking <- data.frame()
                     }
                     rv$district_eval_results <- district_summary
                     
                     # Debug için
                     if(nrow(district_performance) > 0) {
                       message("District performance hesaplandı: ", nrow(district_performance), " satır")
                       rv$district_ranking <- district_performance
                     } else {
                       message("District performance boş!")
                       rv$district_ranking <- NULL
                     } 
                     # Koordinat verilerini sakla (harita için - TGDR olmadan)
                     coord_data <- df %>%
                       select(lon, lat, NAME_1, NAME_2, Ra226_grid, Th232_grid, K40_grid) %>%
                       mutate(
                         point_id = row_number(),
                         district_label = paste(NAME_2, NAME_1, sep = ", ")
                       )
                     
                     # Gerçek değerler varsa ekle
                     if (all(c("Ra226", "Th232", "K40") %in% names(df))) {
                       coord_data <- coord_data %>%
                         mutate(
                           Ra226_obs = df$Ra226,
                           Th232_obs = df$Th232,
                           K40_obs = df$K40
                         )
                     }
                     
                     rv$eval_coordinates <- coord_data
                     
                   } else {
                     # Sadece tahmin sonuçları (gerçek değer yok)
                     
                     # İlçe bazında gruplama (sadece tahminler)
                     district_summary <- df |>
                       group_by(ID_1, NAME_1, ID_2, NAME_2) |>
                       summarise(
                         # Grid tahmin ortalamaları
                         Ra226_pixel_mean = round(mean(Ra226_grid, na.rm = TRUE),2),
                         Th232_pixel_mean = round(mean(Th232_grid, na.rm = TRUE),2),
                         K40_pixel_mean  = round(mean(K40_grid, na.rm = TRUE),2),
                         
                         # Nokta sayısı
                         point_count = n(),
                         .groups = 'drop'
                       )
                     
                     # İlçe tahmin ortalamalarını (genel harita ortalamaları) ekle
                     district_summary <- district_summary |>
                       left_join(ilce_means |> select(ID_2, Ra_mean, Th_mean, K_mean), by = "ID_2") |>
                       rename(
                         Ra226_map_mean = Ra_mean,
                         Th232_map_mean = Th_mean,
                         K40_map_mean   = K_mean
                       )
                     
                     rv$district_eval_results <- district_summary
                     rv$district_metrics <- NULL
                     rv$district_ranking <- NULL
                     
                     # Koordinat verilerini sakla (harita için)
                     coord_data <- df %>%
                       select(lon, lat, NAME_1, NAME_2, Ra226_grid, Th232_grid, K40_grid) %>%
                       mutate(
                         point_id = row_number(),
                         district_label = paste(NAME_2, NAME_1, sep = ", ")
                       )
                     
                     rv$eval_coordinates <- coord_data
                   }
                   
                   incProgress(1.0, detail = if (lang()=="tr") "Tamamlandı!" else "Completed!")
                   showNotification(trn(lang(),"evaluation_complete"), type="default")
                 })
  })
  
  # ---- OUTPUTS ----
  # Metrik tabloları
  output$district_metrics_table <- renderTable({
    req(rv$district_metrics)
    rv$district_metrics
  }, digits = 3)
  
  output$district_metrics_inside_table <- renderTable({
    req(rv$district_metrics_inside)
    rv$district_metrics_inside
  }, digits = 3)
  
  output$point_metrics_table <- renderTable({
    req(rv$point_metrics)
    rv$point_metrics
  }, digits = 3)
  
  output$point_metrics_inside_table <- renderTable({
    req(rv$point_metrics_inside)
    rv$point_metrics_inside
  }, digits = 3)
  
  output$tgdr_metrics_table <- renderTable({
    req(rv$tgdr_metrics)
    rv$tgdr_metrics
  }, digits = 3)
  
  output$tgdr_metrics_inside_table <- renderTable({
    req(rv$tgdr_metrics_inside)
    rv$tgdr_metrics_inside
  }, digits = 3)
  
  # YENİ EKLENECEK OLANLAR
  output$district_stat_tests_inside <- renderTable({
    req(rv$district_stat_tests_inside)
    rv$district_stat_tests_inside
  })
  
  output$point_stat_tests_inside <- renderTable({
    req(rv$point_stat_tests_inside)
    rv$point_stat_tests_inside
  })
  
  output$tgdr_stat_tests_inside <- renderTable({
    req(rv$tgdr_stat_tests_inside)
    rv$tgdr_stat_tests_inside
  })

  # İlçe performans sıralaması tablosu
  output$district_ranking_table <- DT::renderDataTable({
    req(rv$district_ranking)
    
    # Debug kontrolü
    if(is.null(rv$district_ranking) || nrow(rv$district_ranking) == 0) {
      message("District ranking tablosu boş!")
      return(NULL)
    }
    
    message("District ranking tablosu render ediliyor: ", nrow(rv$district_ranking), " satır")
    
    # Sütun isimlerini dile göre ayarla
    ranking_data <- rv$district_ranking
    if (lang() == "tr") {
      colnames(ranking_data) <- c("Sıra", "İl", "İlçe", "Nokta Sayısı", "Ortalama MAPE", 
                                  "Kategori", "Ra226 MAPE", "Th232 MAPE", "K40 MAPE")
    } else {
      colnames(ranking_data) <- c("Rank", "Province", "District", "Point Count", "Mean MAPE",
                                  "Category", "Ra226 MAPE", "Th232 MAPE", "K40 MAPE")
    }
    
    DT::datatable(ranking_data, 
                  options = list(pageLength = 15, scrollX = TRUE),
                  selection = 'single') %>%
      DT::formatRound(columns = c(if(lang()=="tr") 'Ortalama MAPE' else 'Mean MAPE',
                                  'Ra226 MAPE', 'Th232 MAPE', 'K40 MAPE'), digits = 2)
  })
  
  # İstatistiksel testler
  output$district_stat_tests <- renderTable({
    req(rv$district_stat_tests)
    rv$district_stat_tests
  })
  
  output$point_stat_tests <- renderTable({
    req(rv$point_stat_tests)
    rv$point_stat_tests
  })
  
  output$tgdr_stat_tests <- renderTable({
    req(rv$tgdr_stat_tests)
    rv$tgdr_stat_tests
  })
  
  # Elips istatistikleri
  output$ellipse_stats_table <- renderTable({
    req(rv$district_ellipse_stats)
    rv$district_ellipse_stats
  })
  
  output$point_ellipse_stats_table <- renderTable({
    req(rv$point_ellipse_stats)
    rv$point_ellipse_stats
  })
  
  output$tgdr_ellipse_stats_table <- renderTable({
    req(rv$tgdr_ellipse_stats)
    rv$tgdr_ellipse_stats
  })
  
  # Çapraz doğrulama grafikleri
  output$district_validation_plots <- renderPlot({
    req(rv$district_validation_plots)
    rv$district_validation_plots()
  })
  
  output$point_validation_plots <- renderPlot({
    req(rv$point_validation_plots)
    rv$point_validation_plots()
  })
  
  output$tgdr_validation_plots <- renderPlot({
    req(rv$tgdr_validation_plots)
    rv$tgdr_validation_plots()
  })
  
  # Ana değerlendirme tablosu
  output$district_eval_table <- DT::renderDataTable({
    req(rv$district_eval_results)
    DT::datatable(rv$district_eval_results, options = list(pageLength = 15, scrollX = TRUE), selection = 'multiple')
  })
  
  # Koordinat haritası - TGDR temizlenmiş
  output$coord_eval_map <- renderLeaflet({
    req(rv$eval_coordinates)
    
    coord_data <- rv$eval_coordinates
    
    leaflet(coord_data) %>%
      addProviderTiles(providers$CartoDB.Positron) %>%
      addPolygons(data = ilce_sf, color = "gray", weight = 1, fillOpacity = 0.1, opacity = 0.5) %>%
      addCircleMarkers(
        lng = ~lon, lat = ~lat,
        radius = 6, 
        color = "#3498db",
        fillColor = "#3498db", 
        fillOpacity = 0.7,
        stroke = TRUE,
        weight = 2,
        popup = ~paste0(
          "<b>", district_label, "</b><br>",
          "Koordinat: ", round(lon, 4), ", ", round(lat, 4), "<br>",
          "<hr><b>", if(lang()=="tr") "Tahmin Değerleri:" else "Predicted Values:", "</b><br>",
          "Ra-226: ", round(Ra226_grid, 2), " Bq/kg<br>",
          "Th-232: ", round(Th232_grid, 2), " Bq/kg<br>",
          "K-40: ", round(K40_grid, 2), " Bq/kg"
        ),
        layerId = ~point_id
      ) %>%
      addScaleBar(position = "bottomleft",
                  options = list(metric = TRUE, imperial = FALSE, maxWidth = 200, updateWhenIdle = TRUE)) %>%
      fitBounds(min(coord_data$lon), min(coord_data$lat), 
                max(coord_data$lon), max(coord_data$lat))
  })
  
  # Tablo seçiminde haritayı güncelle - TGDR temizlenmiş
  observeEvent(input$district_eval_table_rows_selected, {
    req(rv$eval_coordinates)
    req(input$district_eval_table_rows_selected)
    
    coord_data <- rv$eval_coordinates
    selected_districts <- rv$district_eval_results[input$district_eval_table_rows_selected, ]
    
    # Seçili ilçelerdeki koordinatları belirle
    selected_coords <- coord_data %>%
      filter(NAME_2 %in% selected_districts$NAME_2)
    
    non_selected_coords <- coord_data %>%
      filter(!NAME_2 %in% selected_districts$NAME_2)
    
    # Haritayı güncelle
    leafletProxy("coord_eval_map") %>%
      clearMarkers() %>%
      # Seçilmeyen koordinatlar (gri)
      addCircleMarkers(
        data = non_selected_coords,
        lng = ~lon, lat = ~lat,
        radius = 5, 
        color = "#95a5a6",
        fillColor = "#95a5a6", 
        fillOpacity = 0.5,
        stroke = TRUE,
        weight = 1,
        popup = ~paste0(
          "<b>", district_label, "</b><br>",
          "Koordinat: ", round(lon, 4), ", ", round(lat, 4), "<br>",
          "<hr><b>", if(lang()=="tr") "Tahmin Değerleri:" else "Predicted Values:", "</b><br>",
          "Ra-226: ", round(Ra226_grid, 2), " Bq/kg<br>",
          "Th-232: ", round(Th232_grid, 2), " Bq/kg<br>",
          "K-40: ", round(K40_grid, 2), " Bq/kg"
        ),
        layerId = ~point_id
      ) %>%
      # Seçili koordinatlar (kırmızı)
      addCircleMarkers(
        data = selected_coords,
        lng = ~lon, lat = ~lat,
        radius = 8, 
        color = "#e74c3c",
        fillColor = "#e74c3c", 
        fillOpacity = 0.9,
        stroke = TRUE,
        weight = 3,
        popup = ~paste0(
          "<b>", district_label, " (", if(lang()=="tr") "SEÇİLİ" else "SELECTED", ")</b><br>",
          "Koordinat: ", round(lon, 4), ", ", round(lat, 4), "<br>",
          "<hr><b>", if(lang()=="tr") "Tahmin Değerleri:" else "Predicted Values:", "</b><br>",
          "Ra-226: ", round(Ra226_grid, 2), " Bq/kg<br>",
          "Th-232: ", round(Th232_grid, 2), " Bq/kg<br>",
          "K-40: ", round(K40_grid, 2), " Bq/kg"
        ),
        layerId = ~point_id
      )
  })
  
  # Görünürlük kontrolleri
  output$show_district_metrics <- reactive({ 
    !is.null(rv$district_metrics) 
  })
  outputOptions(output, "show_district_metrics", suspendWhenHidden = FALSE)
  
  output$show_point_metrics <- reactive({ 
    !is.null(rv$point_metrics) 
  })
  outputOptions(output, "show_point_metrics", suspendWhenHidden = FALSE)
  
  output$show_tgdr_metrics <- reactive({ 
    !is.null(rv$tgdr_metrics) 
  })
  outputOptions(output, "show_tgdr_metrics", suspendWhenHidden = FALSE)
  
  output$show_stat_tests <- reactive({ # !is.null(rv$district_stat_tests) && !is.null(rv$point_stat_tests) && !is.null(rv$tgdr_stat_tests)
   FALSE 
  })
  outputOptions(output, "show_stat_tests", suspendWhenHidden = FALSE)
  
  output$show_district_plots <- reactive({ 
    !is.null(rv$district_validation_plots) 
  })
  outputOptions(output, "show_district_plots", suspendWhenHidden = FALSE)
  
  output$show_point_plots <- reactive({ 
    !is.null(rv$point_validation_plots) 
  })
  outputOptions(output, "show_point_plots", suspendWhenHidden = FALSE)
  
  output$show_tgdr_plots <- reactive({ 
    !is.null(rv$tgdr_validation_plots) 
  })
  outputOptions(output, "show_tgdr_plots", suspendWhenHidden = FALSE)
  
  output$show_ellipse_stats <- reactive({ 
    !is.null(rv$district_ellipse_stats) 
  })
  outputOptions(output, "show_ellipse_stats", suspendWhenHidden = FALSE)
  
  output$show_district_ranking <- reactive({ 
    result <- !is.null(rv$district_ranking) && nrow(rv$district_ranking) > 0
    message("District ranking görünürlük: ", result)
    return(result)
  })
  outputOptions(output, "show_district_ranking", suspendWhenHidden = FALSE)
  
  output$show_coord_map <- reactive({ 
    !is.null(rv$eval_coordinates) && nrow(rv$eval_coordinates) > 0
  })
  outputOptions(output, "show_coord_map", suspendWhenHidden = FALSE)
  
  # Download handlers
  output$download_eval_results <- downloadHandler(
    filename = function() {
      paste0(if (lang()=="tr") "kapsamli_ilce_degerlendirme_" else "comprehensive_district_evaluation_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      wb <- openxlsx::createWorkbook()
      
      # Excel sayfaları oluştur
      if (!is.null(rv$district_eval_results)) {
        openxlsx::addWorksheet(wb, if(lang()=="tr") "İlçe Sonuçları" else "District Results")
        openxlsx::writeData(wb, if(lang()=="tr") "İlçe Sonuçları" else "District Results", rv$district_eval_results)
      }
      
      if (!is.null(rv$district_ranking)) {
        openxlsx::addWorksheet(wb, if(lang()=="tr") "Performans Sıralaması" else "Performance Ranking")
        openxlsx::writeData(wb, if(lang()=="tr") "Performans Sıralaması" else "Performance Ranking", rv$district_ranking)
      }
      
      if (!is.null(rv$eval_coordinates)) {
        openxlsx::addWorksheet(wb, if(lang()=="tr") "Koordinatlar" else "Coordinates")
        openxlsx::writeData(wb, if(lang()=="tr") "Koordinatlar" else "Coordinates", rv$eval_coordinates)
      }
      
      if (!is.null(rv$district_metrics)) {
        openxlsx::addWorksheet(wb, if(lang()=="tr") "İlçe Metrikleri" else "District Metrics")
        openxlsx::writeData(wb, if(lang()=="tr") "İlçe Metrikleri" else "District Metrics", rv$district_metrics, startRow = 1)
        
        if (!is.null(rv$district_metrics_inside)) {
          openxlsx::writeData(wb, if(lang()=="tr") "İlçe Metrikleri" else "District Metrics", 
                              data.frame(Section = if(lang()=="tr") "Elips İçi Metrikler (95%)" else "Inside Ellipse Metrics (95%)"), 
                              startRow = nrow(rv$district_metrics) + 3)
          openxlsx::writeData(wb, if(lang()=="tr") "İlçe Metrikleri" else "District Metrics", 
                              rv$district_metrics_inside, startRow = nrow(rv$district_metrics) + 4)
        }
      }
      
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  output$download_eval_plots <- downloadHandler(
    filename = function() {
      paste0(if (lang()=="tr") "capraz_dogrulama_grafikleri_" else "cross_validation_plots_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      pdf(file, width = 16, height = 12)
      
      if (!is.null(rv$district_validation_plots)) {
        rv$district_validation_plots()
        title(main = if(lang()=="tr") "İlçe Bazlı Çapraz Doğrulama" else "District-Based Cross Validation", 
              line = -2, outer = TRUE, cex.main = 1.5)
      }
      
      if (!is.null(rv$point_validation_plots)) {
        plot.new()
        rv$point_validation_plots()
        title(main = if(lang()=="tr") "Nokta Bazlı Çapraz Doğrulama" else "Point-Based Cross Validation", 
              line = -2, outer = TRUE, cex.main = 1.5)
      }
      
      if (!is.null(rv$tgdr_validation_plots)) {
        plot.new()
        rv$tgdr_validation_plots()
        title(main = if(lang()=="tr") "TGDR Çapraz Doğrulama" else "TGDR Cross Validation", 
              line = -2, outer = TRUE, cex.main = 1.5)
      }
      
      dev.off()
    }
  )
  
  output$download_eval_plots_svg <- downloadHandler(
    filename = function() {
      paste0(if (lang()=="tr") "capraz_dogrulama_svg_" else "cross_validation_svg_", Sys.Date(), ".zip")
    },
    content = function(file) {
      # Geçici dizin ve SVG yolları
      tmpdir <- tempfile("svg_export")
      dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
      svg_paths <- character(0)
      
      # 1) İlçe bazlı plot (varsa)
      if (!is.null(rv$district_validation_plots)) {
        f <- file.path(tmpdir, if (lang()=="tr") "ilce_bazli.svg" else "district_based.svg")
        grDevices::svg(filename = f, width = 16, height = 12, onefile = FALSE)
        rv$district_validation_plots()
        dev.off()
        svg_paths <- c(svg_paths, f)
      }
      
      # 2) Nokta bazlı plot (varsa)
      if (!is.null(rv$point_validation_plots)) {
        f <- file.path(tmpdir, if (lang()=="tr") "nokta_bazli.svg" else "point_based.svg")
        grDevices::svg(filename = f, width = 16, height = 12, onefile = FALSE)
        rv$point_validation_plots()
        dev.off()
        svg_paths <- c(svg_paths, f)
      }
      
      # 3) TGDR plot (varsa)
      if (!is.null(rv$tgdr_validation_plots)) {
        f <- file.path(tmpdir, "tgdr.svg")
        grDevices::svg(filename = f, width = 12, height = 9, onefile = FALSE)
        rv$tgdr_validation_plots()
        dev.off()
        svg_paths <- c(svg_paths, f)
      }
      
      # Dosyaları ZIP'e koy
      owd <- setwd(tmpdir)
      on.exit(setwd(owd), add = TRUE)
      
      # 'zip' paketi varsa onu kullan, yoksa utils::zip
      if (requireNamespace("zip", quietly = TRUE)) {
        zip::zipr(zipfile = file, files = basename(svg_paths))
      } else {
        utils::zip(zipfile = file, files = basename(svg_paths))
      }
    }
  )
  
}
#options(shiny.trace = TRUE)
# -------------------- APP ----------------------------------------------------
shinyApp(ui, server)

#save.image("calisan_shine_8.rdata")



