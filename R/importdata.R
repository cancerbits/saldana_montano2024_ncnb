
###############################################################################################################################################
#                                   DEFINING COLORS
###############################################################################################################################################

allcolors=list(ordclusters = c(`1` = "#fda85f", `2` = "#fcc1e2", `3` = "#ffe16c", 
`4` = "#fdb4ae", `5` = "#bc74ba", `6` = "#b3d266", `7` = "#beaed7", 
`8` = "#ff22ff", `9` = "#bebb26", `10` = "#80b21f", `11` = "#fff3b0", 
`12` = "#998d96", `13` = "#fb80be", `14` = "#ccdfc2", `15` = "#8dd413", 
`16` = "#8dc7c4", `17` = "#998d96", `18` = "#fb746f", `19` = "#d9cdd6", 
`20` = "#80a5d0"), ccondition = c(cWT = "#777777", c17q = "#00AD67", 
c17q1q = "orange", c17q1qMYCN = "magenta"), stage = c(D0 = "#9E0142", 
D3 = "#D53E4F", D6 = "#F46D43", D9 = "#FDAE61", D10 = "#FEE08B", 
D12 = "#ABDDA4", D14 = "#66C2A5", D19 = "#3288BD", D28 = "#5E4FA2"
), orig.ident = c("#fda85f", "#fcc1e2", "#ffe16c", "#fdb4ae", 
"#bc74ba", "#b3d266", "#beaed7", "#ff22ff", "#bebb26", "#80b21f", 
"#fff3b0", "#998d96", "#fb80be", "#ccdfc2", "#8dd413", "#8dc7c4", 
"#998d96", "#fb746f", "#d9cdd6", "#80a5d0"), mordclusters = c(M1 = "#fda85f", 
M2 = "#fcc1e2", M3 = "#ffe16c", M4 = "#fdb4ae", M5 = "#bc74ba", 
M6 = "#b3d266", M7 = "#beaed7", M8 = "#ff22ff", M9 = "#bebb26", 
M10 = "#80b21f", M11 = "#fff3b0", M12 = "#998d96", M13 = "#fb80be", 
M14 = "#ccdfc2", M15 = "#8dd413", M16 = "#8dc7c4", M17 = "#998d96", 
M18 = "#fb746f", M19 = "#d9cdd6", M20 = "#80a5d0"), Phase = c(G1 = "#7fc97f", 
S = "#beaed4", G2M = "#fdc086"), wtrclusters = c(C1 = "#8dd3c7", 
C2 = "#ffffb3", C3 = "#bebada", C4 = "#fb8072", C5 = "#80b1d3", 
C6 = "#fdb462", C7 = "#b3de69", C8 = "#fccde5", C9 = "#d9d9d9", 
C10 = "#bc80bd", C11 = "#ccebc5", C12 = "#ffed6f", C13 = "#FF7BBAFF"
), wtclusters_arranged = c(`1` = "#8dd3c7", `2` = "#ffffb3", 
`3` = "#bebada", `4` = "#fb8072", `5` = "#80b1d3", `6` = "#fdb462", 
`7` = "#b3de69", `8` = "#fccde5", `9` = "#d9d9d9", `10` = "#bc80bd", 
`11` = "#ccebc5", `12` = "#ffed6f", `13` = "#FF7BBAFF"), select = c(C1 = NA, 
C2 = NA, C3 = NA, C4 = NA, C5 = NA, C6 = NA, C7 = "red", C8 = NA, 
C9 = NA, C10 = NA, C11 = NA, C12 = NA, C13 = NA), TF = c(`TRUE` = "red", 
`FALSE` = "grey"), ceurat_clusters4 = c(C1 = "#bebada", C2 = "#fccde5", 
C3 = "#80b1d3", C4 = "#FEEDC3FF", C5 = "#ccebc5", C6 = "#FF7BBAFF", 
C7 = "#E6F288FF", C8 = "#B5C8FAFF", C9 = "#ffed6f", C10 = "#FDA6D8FF", 
C11 = "#999999", C12 = "#d9d9d9", C13 = "#b3de69", C14 = "#fdb462", 
C15 = "#bc80bd", C16 = "#B4F9FDFF", C17 = "#fb8072", C18 = "#EE8B6EFF", 
C19 = "#8dd3c7", C20 = "#ffffb3", C21 = "#bebada", C22 = "#fccde5", 
C23 = "#80b1d3", C24 = "#FEEDC3FF", C25 = "#ccebc5", C26 = "#FF7BBAFF", 
C27 = "#E6F288FF", C28 = "#B5C8FAFF", C29 = "#ffed6f", C30 = "#FDA6D8FF", 
C31 = "#999999", C32 = "#d9d9d9", C33 = "#b3de69", C34 = "#fdb462", 
C35 = "#bc80bd", C36 = "#B4F9FDFF", C37 = "#fb8072", C38 = "#EE8B6EFF", 
C39 = "#8dd3c7", C40 = "#ffffb3", C41 = "#bebada", C42 = "#fccde5", 
C43 = "#80b1d3", C44 = "#FEEDC3FF", C45 = "#ccebc5", C46 = "#FF7BBAFF", 
C47 = "#E6F288FF", C48 = "#B5C8FAFF", C49 = "#ffed6f"), WT.clusters = c(C1 = "#03C78F", 
C2 = "#86D36C", C3 = "#CB5CBB", C4 = "#F0C4E5", C5 = "#086625", 
C6 = "#19FAFD", C7 = "#9FE6B2", C8 = "#C190AA", C9 = "#9E9C79", 
C10 = "#F4BE67", C11 = "#4B5702", C12 = "#92A0CD", C13 = "#55B384", 
C14 = "#CD0252"), WT.clusters2 = c(`1` = "#8dd3c7", `2` = "#ffffb3", 
`3` = "#bebada", `4` = "#fb8072", `5` = "#80b1d3", `6` = "#fdb462", 
`7` = "#b3de69", `8` = "#fccde5", `9` = "#d9d9d9", `10` = "#bc80bd", 
`11` = "#ccebc5", `12` = "#ffed6f"), WT.clusters3 = c(`1` = "#8dd3c7", 
`2` = "#ffffb3", `3` = "#bebada", `4` = "#fb8072", `5` = "#80b1d3", 
`6` = "#fdb462", `7` = "#b3de69", `8` = "#fccde5", `9` = "#d9d9d9", 
`10` = "#bc80bd", `11` = "#ccebc5"), WT.clusters4 = c(`1` = "#bebada", 
`2` = "#fccde5", `3` = "#80b1d3", `4` = "#FEEDC3FF", `5` = "#ccebc5", 
`6` = "#FF7BBAFF", `7` = "#E6F288FF", `8` = "#B5C8FAFF", `9` = "#ffed6f", 
`10` = "#FDA6D8FF", `11` = "#999999"), WT.clusters5 = c(`1` = "#bebada", 
`2` = "#fccde5", `3` = "#80b1d3", `4` = "#FEEDC3FF", `5` = "#ccebc5", 
`6` = "#FF7BBAFF", `7` = "#E6F288FF", `8` = "#B5C8FAFF", `9` = "#ffed6f"
), WT.clusters6 = c(`1` = "#bebada", `2` = "#fccde5", `3` = "#80b1d3", 
`4` = "#FEEDC3FF", `5` = "#ccebc5", `6` = "#FF7BBAFF", `7` = "#E6F288FF", 
`8` = "#B5C8FAFF", `9` = "#ffed6f", `10` = "#FDA6D8FF"), WT.clusters7 = c(`1` = "#bebada", 
`2` = "#fccde5", `3` = "#80b1d3", `4` = "#FEEDC3FF", `5` = "#ccebc5", 
`6` = "#FF7BBAFF", `7` = "#E6F288FF", `8` = "#B5C8FAFF", `9` = "#ffed6f", 
`10` = "#FDA6D8FF", `11` = "#999999"), type = list(SCP = "#FDA6D8FF", 
    mesenchyme = "#8dd3c7", sympathoblasts = "#fdb462", other = "#999999"), 
    cutoff_clusters_arranged = c(`1` = "#bebada", `2` = "#fccde5", 
    `3` = "#80b1d3", `4` = "#FEEDC3FF", `5` = "#ccebc5", `6` = "#FF7BBAFF", 
    `7` = "#E6F288FF", `8` = "#B5C8FAFF", `9` = "#ffed6f", `10` = "#FDA6D8FF", 
    `11` = "#999999", `12` = "#d9d9d9", `13` = "#b3de69"), cutoff_clusters = c(`1` = "#bebada", 
    `2` = "#fccde5", `3` = "#80b1d3", `4` = "#FEEDC3FF", `5` = "#ccebc5", 
    `6` = "#FF7BBAFF", `7` = "#E6F288FF", `8` = "#B5C8FAFF", 
    `9` = "#ffed6f", `10` = "#FDA6D8FF", `11` = "#999999", `12` = "#d9d9d9", 
    `13` = "#b3de69", `14` = "#fdb462", `15` = "#bc80bd", `16` = "#B4F9FDFF", 
    `17` = "#fb8072", `18` = "#EE8B6EFF"), ceurat_mutant_clusters = c(M1 = "#bebada", 
    M2 = "#fccde5", M3 = "#80b1d3", M4 = "#FEEDC3FF", M5 = "#ccebc5", 
    M6 = "#FF7BBAFF", M7 = "#E6F288FF", M8 = "#B5C8FAFF", M9 = "#ffed6f", 
    M10 = "#FDA6D8FF", M11 = "#999999", M12 = "#d9d9d9", M13 = "#b3de69", 
    M14 = "#fdb462", M15 = "#bc80bd", M16 = "#B4F9FDFF", M17 = "#AAAAAA", 
    M18 = "#EE8B6EFF", M19 = "#8dd3c7", M20 = "#ffffb3", M21 = "#bebada", 
    M22 = "#fccde5", M23 = "#80b1d3", M24 = "#FEEDC3FF", M25 = "#ccebc5"
    ), mapfun_fate2 = c(HSC_and_immune = "#bebada", intermediate_mesoderm = "#fccde5", 
    sympathoblasts = "#fdb462", endothelium = "#FEEDC3FF", cortex = "#ccebc5", 
    kidney = "#E6F288FF", mesenchyme = "#8dd3c7", chromaffin = "#ffed6f", 
    SCP = "#FDA6D8FF", melanocytes = "#999999", liver = "#80b1d3"
    ), mapfun_ordclusters = c(`2` = "#bebada", `4` = "#fccde5", 
    `19` = "#80b1d3", `13` = "#FEEDC3FF", `8` = "#ccebc5", `11` = "#FF7BBAFF", 
    `7` = "#E6F288FF", `1` = "#B5C8FAFF", `20` = "#ffed6f", `14` = "#FDA6D8FF", 
    `18` = "#999999", `15` = "#d9d9d9", `6` = "#b3de69", `12` = "#fdb462", 
    `16` = "#bc80bd", `5` = "#B4F9FDFF", `9` = "#fb8072", `3` = "#EE8B6EFF"
    ), clusters = c(`1` = "#bebada", `2` = "#fccde5", `3` = "#80b1d3", 
    `4` = "#FEEDC3FF", `5` = "#ccebc5", `6` = "#FF7BBAFF", `7` = "#E6F288FF", 
    `8` = "#B5C8FAFF", `9` = "#ffed6f", `10` = "#FDA6D8FF", `11` = "#999999"
    ), WT.clusters8 = c(C1 = "#bebada", C2 = "#fccde5", C3 = "#80b1d3", 
    C4 = "#FEEDC3FF", C5 = "#ccebc5", C6 = "#FF7BBAFF", C7 = "#E6F288FF", 
    C8 = "#B5C8FAFF", C9 = "#ffed6f", C10 = "#FDA6D8FF", C11 = "#999999"
    ), roi = c(R1 = "#BA6ED1", R2 = "#223658", R3 = "#536C98", 
    R4 = "#0EC0B4", R5 = "#1BC242", R6 = "#2648E8", R7 = "#086431", 
    R8 = "#51912F"), replicate = c(R1 = "red", R2 = "blue", R3 = "yellow", 
    R4 = "green", R5 = "purple"), customcols = c(M1 = "grey", 
    M3 = "grey", M5 = "grey", M6 = "grey", M13 = "grey", M18 = "grey", 
    M8 = "grey", M4 = "grey", M16 = "grey", M2 = "grey", M15 = "#ff00ff", 
    M10 = "grey", M7 = "#ba00ff", M14 = "grey", M22 = "grey", 
    M24 = "#e83f47", M12 = "grey", M23 = "grey", M20 = "#8200ff", 
    M19 = "grey", M9 = "#f1885a", M11 = "grey", M17 = "#FDB461", 
    M21 = "grey", M25 = "grey"), WT.main.clusters = c(C1 = "#bebada", 
    C2 = "#fccde5", C3 = "#80b1d3", C4 = "#FEEDC3FF", C5 = "#ccebc5", 
    C6 = "#FF7BBAFF", C7 = "#E6F288FF", C8 = "#B5C8FAFF", C9 = "#ffed6f", 
    C10 = "#FDA6D8FF", C11 = "#999999", C_ = "grey"), megaclusters6 = c(STEMCELLS = "#99999977", 
    SYM = "#FDB46277", MESSYM = "#029ED077", MES = "#8DD3C777", 
    EARLYMES = "steelblue", SCPMES = "#EE365A", SCP = "#FDA6D877", 
    SCPSYM = "black", LATESENSORY = "#74479C77", EARLYSENSORY = "#A282C477", 
    NC = "#08A94B77", HD17q1q = "#F2EA2F77", HDMYCN_D9 = "#FCA1FF", 
    HDMYCN_D14 = "#FF00FF", HDMYCN_D19 = "#B500BA", NMP = "#0C91F477", 
    ENDOTHELIAL = "#A82A2377"), patientID = c(Jansky_NB01 = "#9E9C79", 
    Jansky_NB08 = "#F4BE67", Jansky_NB14 = "#92A0CD", dong_T162 = "#55B384", 
    dong_T200 = "#CD0252", dong_T230 = "#6D22B0", Fetahu_M2 = "#8B1750", 
    Fetahu_M4 = "#16EA6E", Fetahu_M3 = "#38724F", Fetahu_M1 = "#231F0B"
    ), mapfun_seurat_clusters = c(C1 = "#F96300", C2 = "#B45A1E", 
    C3 = "#9D314B", C4 = "#B4B850", C5 = "#C4B820", C6 = "#D58FB0", 
    C7 = "#9149B7", C8 = "#4E2A90", C9 = "#994AD5", C10 = "#369B85", 
    C11 = "#E93162", C12 = "#3E168B", C13 = "#7C8014", C14 = "#9B0AA2"
    ), marker_source = c(D0 = "#d7191c", D3 = "#fdae61", D9 = "#eeeebf", 
    D14 = "#abdda4", D19 = "#2b83ba", clinical = "purple", invitro_mutant = "magenta"
    ), adameykotype = c(HSC_and_immune = "#bebada", intermediate_mesoderm = "#CE9B70", 
    sympathoblasts = "#fdb462", endothelium = "#B64BD4", cortex = "#87b1d2", 
    HSC_and_immune = "firebrick", kidney = "#E6F288FF", mesenchyme = "#8dd3c7", 
    chromaffin = "#438900", SCP = "#FDA6D8FF", melanocytes = "#999999", 
    liver = "#FF7BBAFF"), mapfun_ccondition = c(cWT = "#777777", 
    c17q = "#00AD67", c17q1q = "orange", c17q1qMYCN = "magenta"
    ), null = c(one = "#FFFFFF00"), ccondition_old = c(cWT = "#3288BD", 
    c17q = "#00AD67", c17q1q = "orange", c17q1qMYCN = "magenta"
    ), combinedcat = c(D0_c17q1q_other = "#999999", D0_c17q_other = "#999999", 
    D9_c17q1q_SCP = "#08A94B", D0_cWT_other = "#999999", D9_c17q1qMYCN_other = "#FCA1FF", 
    D9_c17q_SCP = "#08A94B", D9_cWT_SCP = "#08A94B", D19_c17q_sympathoblasts = "#FDB462", 
    D9_c17q1q_other = "#F2EA2F", D19_cWT_mesenchyme = "#8DD3C7", 
    D9_cWT_sympathoblasts = "#74479C", D19_cWT_sympathoblasts = "#FDB462", 
    D9_cWT_other = "#A282C4", D14_c17q_mesenchyme = "#8DD3C7", 
    D3_c17q1q_other = "#0C91F4", D14_cWT_other = "black", D9_c17q1qMYCN_SCP = "#FCA1FF", 
    D19_c17q_mesenchyme = "#8DD3C7", D9_c17q_other = "#A282C4", 
    D19_cWT_other = "#029ED0", D14_c17q_other = "black", D19_c17q1qMYCN_other = "#B500BA", 
    D19_cWT_SCP = "#FDA6D8", D14_c17q1qMYCN_other = "#FF00FF", 
    D19_c17q1q_other = "#EE365A", D14_cWT_mesenchyme = "#8DD3C7", 
    D19_c17q1q_SCP = "#FDA6D8", D19_c17q1q_mesenchyme = "#8DD3C7", 
    D14_cWT_SCP = "#FDA6D8", D14_c17q_SCP = "#FDA6D8", D14_c17q_sympathoblasts = "#FDB462", 
    D14_c17q1q_sympathoblasts = "#FDB462", D9_c17q1q_sympathoblasts = "#FDB462", 
    D19_c17q1qMYCN_SCP = "#FDA6D8", D14_c17q1q_SCP = "black", 
    D14_cWT_sympathoblasts = "#FDB462", D19_c17q_other = "#AF781B", 
    D19_c17q_SCP = "#979E34", D3_c17q_other = "#8F46E4", D19_c17q1qMYCN_sympathoblasts = "#CBDDB5", 
    D3_cWT_other = "#3FE10F", D19_c17q1q_sympathoblasts = "black", 
    D9_c17q_mesenchyme = "#8DD3C7", D19_c17q1qMYCN_mesenchyme = "#8DD3C7", 
    D9_c17q_sympathoblasts = "#4F0439", D9_c17q1q_mesenchyme = "#8DD3C7", 
    D14_c17q1q_mesenchyme = "#8DD3C7", D14_c17q1qMYCN_SCP = "#C17DD3", 
    D14_c17q1q_other = "#974477", D9_cWT_mesenchyme = "#8DD3C7", 
    D14_c17q1qMYCN_sympathoblasts = "#88A8E8", D9_c17q1qMYCN_mesenchyme = "#5D532C", 
    D0_c17q1q_SCP = "#2C8868", D14_c17q1qMYCN_mesenchyme = "#ED4FBE", 
    D9_c17q1qMYCN_sympathoblasts = "#6B3D3C", D0_c17q_mesenchyme = "#E9B70B", 
    D0_c17q_SCP = "#DEBCFB", D0_cWT_mesenchyme = "#91900B", D3_cWT_mesenchyme = "#344587"
    ), wtumap.sctharmony_snn_res.0.4.arranged = c(C1 = "#ADC7B8", 
    C2 = "#E26F6C", C3 = "#8ADC4C", C4 = "#0BAE77", C5 = "#866722", 
    C6 = "#16272F", C7 = "#3CB443", C8 = "#5E8576", C9 = "#84EBE1", 
    C10 = "#6FCC44", C11 = "#2DF9A8", C12 = "#F472A7", C13 = "#31C44C", 
    C14 = "#AD32C3", C15 = "#82D594", C16 = "#CCA2EA", C17 = "#CFF3B5", 
    C18 = "#38C02C", C19 = "#3C08F5", C20 = "#03D281", C21 = "#BC026D", 
    C22 = "#BF8913"), petrolgold = function (n) 
    {
        x <- ramp(seq.int(0, 1, length.out = n))
        if (ncol(x) == 4L) 
            rgb(x[, 1L], x[, 2L], x[, 3L], x[, 4L], maxColorValue = 255)
        else rgb(x[, 1L], x[, 2L], x[, 3L], maxColorValue = 255)
    }, conditiondemux2 = c("#CBDE73", "#445289", "#97032E", "#86B706"
    ), dsname = c(G1 = "#A2FF37", G3 = "#F7DD86", G4 = "#BADAEA", 
    G5 = "#3DDA3D", G6 = "#57321A", G7 = "#665CFD", G8 = "#E92C90", 
    G9 = "#D9D580", G10 = "#9881FF", G11 = "#2F5AC8", G12 = "#8893EC", 
    G13 = "#F040AD", G14 = "#ABB092", G21 = "#7AAD43", G15 = "#DB5BAB", 
    G22 = "#6EA851", G16 = "#D5092F", G23 = "#44ABC2", G17 = "#FB085D", 
    G24 = "#DEEFEA", G18 = "#27CFF9", G25 = "#3256FB", G19 = "#9BC62E", 
    G26 = "#12EA98", G20 = "#1C18BC", G27 = "#038C5A"), umap.fullscvi8_nngraph_res.0.2.arranged = c(M1 = "#C569AB", 
    M2 = "#7CE392", M3 = "#F63C0C", M4 = "#A2AD60", M5 = "#4AFD48", 
    M6 = "#436FCF", M7 = "#D3049F", M8 = "#D70749", M9 = "#A6A4CA", 
    M10 = "#725557", M11 = "#BA6832", M12 = "#E8E660", M13 = "#7E3A2A", 
    M14 = "#43100E", M15 = "#F35034", M16 = "#8D4FD2", M17 = "#760CF5", 
    M18 = "#8E75A2"), mutant.clusters = c(M1 = "#C9F04D", M2 = "#A85748", 
    M3 = "#70E3D9", M4 = "#46359E", M5 = "#9D5D7C", M6 = "#891F8D", 
    M7 = "#8976F2", M8 = "#83DA92", M9 = "#9133A6", M10 = "#B6A126", 
    M11 = "#B33702", M12 = "#692D5C", M13 = "#DD4085", M14 = "#D9AC6D", 
    M15 = "#AD0869", M16 = "#8722DD"), mutant.clusters.numbers = c(`1` = "#C9F04D", 
    `2` = "#A85748", `3` = "#70E3D9", `4` = "#46359E", `5` = "#9D5D7C", 
    `6` = "#891F8D", `7` = "#8976F2", `8` = "#83DA92", `9` = "#9133A6", 
    `10` = "#B6A126", `11` = "#B33702", `12` = "#692D5C", `13` = "#DD4085", 
    `14` = "#D9AC6D", `15` = "#AD0869", `16` = "#8722DD"), type2 = c(SCP = "#fb8072", 
    mesenchyme = "#8dd3c7", sympathoblasts = "#fdb462", other = "#d9d9d9"
    ), wtumap.fullscvi822_nngraph_res.0.6 = c(`29` = "#E311B0"), 
    mapfun_wtumap.fullscvi822_nngraph_res.0.6 = c(`29` = "#E311B0"), 
    predicted.id = c(C1 = "#198342", C2 = "#993218", C3 = "#233F79", 
    C4 = "#7F92D2", C5 = "#4E7C8E", C6 = "#FD6B01", C7 = "#52B42A", 
    C8 = "#5C62B2", C9 = "#69C277", C10 = "#67C221", C11 = "#767819", 
    C12 = "#F2EBB0", C13 = "#5231A1"), wtumap.fullscvi822_snngraph_res.0.08.arranged = c(C1 = "#AEA936", 
    C2 = "#A98D5A", C3 = "#A7AAF0", C4 = "#32EDE1", C5 = "#3A2FF9", 
    C6 = "#B5CBEE", C7 = "#E7F362", C8 = "#7EB6AD", C9 = "#D9F19F", 
    C10 = "#374566", C11 = "#53BD50", C12 = "#AD1524", C13 = "#C3FE8F"
    ), mapfun_WT.clusters = c(C1 = "#03C78F", C2 = "#86D36C", 
    C3 = "#CB5CBB", C4 = "#F0C4E5", C5 = "#086625", C6 = "#19FAFD", 
    C7 = "#9FE6B2", C8 = "#C190AA", C9 = "#9E9C79", C10 = "#F4BE67", 
    C11 = "#4B5702", C12 = "#92A0CD", C13 = "#55B384", C14 = "#CD0252"
    ), mapped.WT.clusters = c(C1 = "#198342", C2 = "#993218", 
    C3 = "#233F79", C4 = "#7F92D2", C5 = "#4E7C8E", C6 = "#FD6B01", 
    C7 = "#52B42A", C8 = "#5C62B2", C9 = "#69C277", C10 = "#67C221", 
    C11 = "#767819", C12 = "#F2EBB0", C13 = "#5231A1"), WT.clusters.arranged = c(C1 = "#C65733", 
    C2 = "#F57B4D", C3 = "#AFD7C4", C4 = "#485793", C5 = "#89A2B9", 
    C6 = "#330949", C7 = "#EF5FDE", C8 = "#BBF875", C9 = "#AE7369", 
    C10 = "#56397E", C11 = "#F9423A", C12 = "#836BFC", C13 = "#2971E3", 
    C14 = "#36C7AC"), megaclusters = c(STEMCELLS = "#97C964", 
    NMP = "#8E4064", NC = "#AE5D72", NC_UNDIF = "#2496D8", MYCN_EARLY = "#E689ED", 
    SENSORY_CNA = "#0EA8E4", MYCN_MID = "#1E4994", SENSORY = "#054B7B", 
    MES_UNDIF = "#7D8461", SCP = "#760C54", MYCN_SCP = "#F62CD5", 
    MYCN_LATE = "#0FE856", SYM = "#8FC572", MES_SYM = "#20698A", 
    MES = "#C6CB9F", MES_CNA = "#08BE1B"), mapfun_megaclusters = c(STEMCELLS = "#9CE4AA", 
    NMP = "#C9A612", NC = "#81F5B6", NC_UNDIF = "#D04AA8", MYCN_EARLY = "#D1609E", 
    SENSORY_CNA = "#710497", MYCN_MID = "#39DBDB", SENSORY = "#E85D40", 
    MES_UNDIF = "#3C6200", SCP = "#23661B", MYCN_SCP = "#B61A03", 
    MYCN_LATE = "#FAED07", SYM = "#7BA829", MES_SYM = "#4AC0DD", 
    MES = "#B0AB3F", MES_CNA = "#30AD90"), seurat_clusters = c(C1 = "#F53464", 
    C2 = "#70930A", C3 = "#D7575D", C4 = "#066D68", C5 = "#1EC3BC", 
    C6 = "#40BE80", C7 = "#92444B", C8 = "#2FED43", C9 = "#6FEE0C", 
    C10 = "#284BE9", C11 = "#1F8B2A", C12 = "#6B1E73", C13 = "#B89720", 
    C14 = "#498E00"), kamenevatype = c(HSC_and_immune = "#bebada", 
    intermediate_mesoderm = "#fccde5", sympathoblasts = "#fdb462", 
    endothelium = "#FEEDC3FF", cortex = "#ccebc5", kidney = "#E6F288FF", 
    mesenchyme = "#8dd3c7", chromaffin = "#ffed6f", SCP = "#FDA6D8FF", 
    melanocytes = "#999999", liver = "#80b1d3"), mapfun_kamenevatype = c(HSC_and_immune = "#bebada", 
    intermediate_mesoderm = "#fccde5", sympathoblasts = "#fdb462", 
    endothelium = "#FEEDC3FF", cortex = "#ccebc5", kidney = "#E6F288FF", 
    mesenchyme = "#8dd3c7", chromaffin = "#ffed6f", SCP = "#FDA6D8FF", 
    melanocytes = "#999999", liver = "#80b1d3"), mapfun_mutant.clusters = c(M1 = "#C569AB", 
    M2 = "#7CE392", M3 = "#F63C0C", M4 = "#A2AD60", M5 = "#4AFD48", 
    M6 = "#436FCF", M7 = "#D3049F", M8 = "#D70749", M9 = "#A6A4CA", 
    M10 = "#725557", M11 = "#BA6832", M12 = "#E8E660", M13 = "#7E3A2A", 
    M14 = "#43100E", M15 = "#F35034", M16 = "#8D4FD2", M17 = "#760CF5", 
    M18 = "#8E75A2"), wtumap.fullscvi822_nngraph_res.0.4.arranged = c(C1 = "#C32D7A", 
    C2 = "#8201EE", C3 = "#B70A0A", C4 = "#2D2971", C5 = "#CBA2B7", 
    C6 = "#F4760D", C7 = "#E23B13", C8 = "#F83BA7", C9 = "#E04364", 
    C10 = "#D842C2", C11 = "#F8B675", C12 = "#F95A80", C13 = "#FBF3DE", 
    C14 = "#EF9663"), type3 = c(SCP = "#fb8072", mesenchyme = "#8dd3c7", 
    sympathoblasts = "#fdb462", chromaffin = "#C045FF", other = "#d9d9d9"
    ), stageroma = c(D0 = "#7E1900", D3 = "#A0621C", D6 = "#C0A439", 
    D9 = "#E3E086", D10 = "#D1ECC8", D12 = "#79D4D9", D14 = "#479AC5", 
    D19 = "#3165AD", D28 = "#1A3399"), stageoriginal = c(D0 = "black", 
    D3 = "brown", D6 = "darkslateblue", D9 = "blue", D10 = "purple", 
    D12 = "red", D14 = "salmon", D19 = "orange", D28 = "yellow"
    ), stagespectral = c(D0 = "#9E0142", D3 = "#D53E4F", D6 = "#F46D43", 
    D9 = "#FDAE61", D10 = "#FEE08B", D12 = "#ABDDA4", D14 = "#66C2A5", 
    D19 = "#3288BD", D28 = "#5E4FA2"), primary_diagnosis = c(`Neuroblastoma, NOS` = "#EED6D5", 
    Ganglioneuroblastoma = "#AE0E97"), tissue_or_organ_of_origin = c(Unknown = "#23C5FC", 
    `Adrenal gland, NOS` = "#896C0B", `Medulla of adrenal gland` = "#369889", 
    `Posterior mediastinum` = "#09B424", `Unknown primary site` = "#00E716", 
    `Pelvis, NOS` = "#6F178D", `Abdomen, NOS` = "#5B5501", `Peripheral nerves and autonomic nervous system of trunk, NOS` = "#040253", 
    Retroperitoneum = "#4C3326", `Kidney, NOS` = "#8FF728", `Peripheral nerves and autonomic nervous system of thorax` = "#0F0188", 
    `Thorax, NOS` = "#7A1AAD", `Mediastinum, NOS` = "#0BF216", 
    `Pelvic lymph nodes` = "#930AE6", `Overlapping lesion of retroperitoneum and peritoneum` = "#6B75E1", 
    Liver = "#D4E4C4", `Skin of lower limb and hip` = "#756871", 
    `Head, face or neck, NOS` = "#42E19C", `Aortic body and other paraganglia` = "#8F3FB5", 
    `Bone marrow` = "#ABF99E", `Intra-abdominal lymph nodes` = "#385852", 
    `Cortex of adrenal gland` = "#5B9A91", `Other ill-defined sites` = "#20399E", 
    `Connective, subcutaneous and other soft tissues of abdomen` = "#F73CDF"
    ), mapfun_janskytype = c(`Mesenchymal cells` = "#C9DCEB", 
    `Adrenal cortex` = "#8614E6", `late SCPs` = "#0B622E", `connecting Chromaffin cells` = "#0D8456", 
    `late Neuroblasts` = "#5F7A06", `cycling Neuroblasts` = "#3714D1", 
    SCPs = "#BFAC18", Neuroblasts = "#DEB661", `cycling SCPs` = "#FB3E01", 
    Bridge = "#21119B", `Muscle progenitor cells` = "#514ACE", 
    Myofibroblasts = "#B80F85", `Endothelial cells` = "#A4A803", 
    `late Chromaffin cells` = "#C211A0"), umap.fullscvi8_nngraph_res.0.6.arranged = c(M1 = "#C9F04D", 
    M2 = "#A85748", M3 = "#70E3D9", M4 = "#46359E", M5 = "#9D5D7C", 
    M6 = "#891F8D", M7 = "#8976F2", M8 = "#83DA92", M9 = "#9133A6", 
    M10 = "#B6A126", M11 = "#B33702", M12 = "#692D5C", M13 = "#DD4085", 
    M14 = "#D9AC6D", M15 = "#AD0869", M16 = "#8722DD"), chromosome = c(chr1 = "orange", 
    chr2 = "#939393", chr3 = "#D8D8D8", chr4 = "#515051", chr5 = "#939393", 
    chr6 = "#D8D8D8", chr7 = "#515051", chr8 = "#939393", chr9 = "#D8D8D8", 
    chr10 = "#515051", chr11 = "#939393", chr12 = "#D8D8D8", 
    chr13 = "#515051", chr14 = "#939393", chr15 = "#D8D8D8", 
    chr16 = "#515051", chr17 = "#00AD67", chr18 = "#D8D8D8", 
    chr19 = "#515051", chr20 = "#939393", chr21 = "#D8D8D8", 
    chr22 = "#515051"), MYCN_status = c(Amplified = "#FF000099", 
    `Not Amplified` = "#0000FF99"), inss_stage = c(`Stage 1` = "#fbb4ae", 
    `Stage 2-3` = "#b3cde3", `Stage 4` = "#ccebc5", `Stage 4s` = "#decbe4"
    ), study = c(Fetahu2023 = "#8da0cb", Dong2020 = "#e78ac3", 
    Jansky2021 = "#a6d854"))


################################################################################
# other-origin "small" Files found on the input_data folder
################################################################################
coredatapaths=list(
  cnv=list(
    "17q1q.cnr",
    "17q.cnr",
    "H7S14.cnr"
  ),
  gtf="gencode.v40.annotation.gtf.gz", 
  fishtab="Copy of tumor size footprint area micrometer.csv",
  ptmetadata="ptdatasetsmeta_noqc.rds",
  ptcounts="ptdatasetscounts_noqc.rds"
)

fulldatapaths=lapply(coredatapaths, function(x) paste0(datapath, x))

countdataroots=list(
  new="~/mnt_out/processed_data_new/", # (new dataset) final destination of processing with cellranger 7.1.0 run by us. 
  old="~/mnt_out/processed_data/" # (old dataset) final destination of processing with cellranger 7.1.0 july 2023
)

countdatafolders=list(
   new="processed_data_new", # (new dataset) final destination of processing with cellranger 7.1.0 run by us. 
   old="processed_data" # (old dataset) final destination of processing with cellranger 7.1.0 july 2023
)

expbatches=list(
G1 ="old",
G2 ="old",
G3 ="old",
G4 ="old",
G5 ="old",
G6 ="old",
G7 ="old",
G8 ="old",
G9 ="old",
G10="old",
G11="old",
G12="old",
G13="old",
G14="new",
G21="new",
G15="new",
G22="new",
G16="new",
G23="new",
G17="new",
G24="new",
G18="new",
G25="new",
G19="new",
G26="new",
G20="new",
G27="new"
)

################################################################################
# Count matrix files
################################################################################


#updated to the counts obtained from cellranger 7.1.0 20230727
countdatamatrices=list(
G14= "day3rep1/outs/multi/count/raw_feature_bc_matrix.h5",
G21= "day3rep2/outs/multi/count/raw_feature_bc_matrix.h5",
G15= "day6rep1/outs/multi/count/raw_feature_bc_matrix.h5",
G22= "day6rep2/outs/multi/count/raw_feature_bc_matrix.h5",
G16= "day10rep1/outs/multi/count/raw_feature_bc_matrix.h5",
G23= "day10rep2/outs/multi/count/raw_feature_bc_matrix.h5",
G17="day12rep1/outs/multi/count/raw_feature_bc_matrix.h5",
G24="day12rep2/outs/multi/count/raw_feature_bc_matrix.h5",
G18="day14rep1/outs/multi/count/raw_feature_bc_matrix.h5",
G25="day14rep2/outs/multi/count/raw_feature_bc_matrix.h5",
G19="day19rep1/outs/multi/count/raw_feature_bc_matrix.h5",
G26="day19rep2/outs/multi/count/raw_feature_bc_matrix.h5",
G20="day28rep1/outs/multi/count/raw_feature_bc_matrix.h5",
G27="day28rep2/outs/multi/count/raw_feature_bc_matrix.h5"
)



dsname.dayrep.map=list(
  G1= "G1_GEX",
G2="G2_GEX",
G3="G3_GEX",
G4="G4_GEX",
G5="G5_GEX",
G6="G6_GEX",
G7="G7_GEX",
G8="G8_GEX",
G9="G9_GEX",
G10="G10_GEX",
G11="G11_GEX",
G12="G12_GEX",
G13="G13_GEX",
  G14= "day3rep1",
 G21= "day3rep2",
 G15= "day6rep1",
 G22= "day6rep2",
G16= "day10rep1",
G23= "day10rep2",
 G17="day12rep1",
 G24="day12rep2",
 G18="day14rep1",
 G25="day14rep2",
 G19="day19rep1",
 G26="day19rep2",
 G20="day28rep1",
 G27="day28rep2"
)


dsname.fullday.map=list(

 G14= "day3",
 G21= "day3",
 G15= "day6",
 G22= "day6",
G16= "day10",
G23= "day10",
 G17="day12",
 G24="day12",
 G18="day14",
 G25="day14",
 G19="day19",
 G26="day19",
 G20="day28",
 G27="day28"
)


stages=list(
G1= "D0",
G2= "D3",
G3= "D9",
G4= "D14",
G5= "D19",
G6= "D9",
G7= "D14",
G8= "D0",
G9= "D3",
G10="D9",
G11="D19",
G12="D14",
G13=  "D0",
G14= "D3",
G21= "D3",
G15= "D6",
G22= "D6",
G16= "D10",
G23= "D10",
G17="D12",
G24="D12",
G18="D14",
G25="D14",
G19="D19",
G26="D19",
G20="D28",
G27="D28"
)

#old data paths
countdatapathsold=lapply(countdatamatricesold, function(x) paste0(countdataroots$old, x))
#new data paths                 
countdatapathsnew= lapply(countdatamatrices, function(x) paste0(countdataroots$new, x))

allcountpaths= c(countdatapathsold, countdatapathsnew) #add when cellranger outputs are ready
################################################################################
#HTO metadata files
################################################################################

allcountpaths.rearranged=allcountpaths[paste0("G", 1:27)]

demuxfilenamesold=list(
G1=  "G1_GEX/outs/multi/multiplexing_analysis/assignment_confidence_table.csv",
G2=  "G2_GEX/outs/multi/multiplexing_analysis/assignment_confidence_table.csv",
G3=  "G3_GEX/outs/multi/multiplexing_analysis/assignment_confidence_table.csv",
G4=  "G4_GEX/outs/multi/multiplexing_analysis/assignment_confidence_table.csv",
G5=  "G5_GEX/outs/multi/multiplexing_analysis/assignment_confidence_table.csv",
G6=  "G6_GEX/outs/multi/multiplexing_analysis/assignment_confidence_table.csv",
G7=  "G7_GEX/outs/multi/multiplexing_analysis/assignment_confidence_table.csv",
G8=  "G8_GEX/outs/multi/multiplexing_analysis/assignment_confidence_table.csv",
G9=  "G9_GEX/outs/multi/multiplexing_analysis/assignment_confidence_table.csv",
G10="G10_GEX/outs/multi/multiplexing_analysis/assignment_confidence_table.csv",
G11="G11_GEX/outs/multi/multiplexing_analysis/assignment_confidence_table.csv",
G12="G12_GEX/outs/multi/multiplexing_analysis/assignment_confidence_table.csv",
G13="G13_GEX/outs/multi/multiplexing_analysis/assignment_confidence_table.csv")




demuxfilenamesnew=list(
G14=  "day3rep1/outs/multi/multiplexing_analysis/assignment_confidence_table.csv",
G21=  "day3rep2/outs/multi/multiplexing_analysis/assignment_confidence_table.csv",
G15=  "day6rep1/outs/multi/multiplexing_analysis/assignment_confidence_table.csv",
G22=  "day6rep2/outs/multi/multiplexing_analysis/assignment_confidence_table.csv",
G16= "day10rep1/outs/multi/multiplexing_analysis/assignment_confidence_table.csv",
G23= "day10rep2/outs/multi/multiplexing_analysis/assignment_confidence_table.csv",
G17= "day12rep1/outs/multi/multiplexing_analysis/assignment_confidence_table.csv",
G24= "day12rep2/outs/multi/multiplexing_analysis/assignment_confidence_table.csv",
G18= "day14rep1/outs/multi/multiplexing_analysis/assignment_confidence_table.csv",
G25= "day14rep2/outs/multi/multiplexing_analysis/assignment_confidence_table.csv",
G19= "day19rep1/outs/multi/multiplexing_analysis/assignment_confidence_table.csv",
G26= "day19rep2/outs/multi/multiplexing_analysis/assignment_confidence_table.csv",
G20= "day28rep1/outs/multi/multiplexing_analysis/assignment_confidence_table.csv",
G27= "day28rep2/outs/multi/multiplexing_analysis/assignment_confidence_table.csv")





demuxfilenames=list(
G14= "day3rep1_transcriptome/HTO_demux.csv",
G21= "day3rep2_transcriptome/HTO_demux.csv",
G15= "day6rep1_transcriptome/HTO_demux.csv",
G22= "day6rep2_transcriptome/HTO_demux.csv",
G16= "day10rep1_transcriptome/HTO_demux.csv",
G23= "day10rep2_transcriptome/HTO_demux.csv",
G17="day12rep1_transcriptome/HTO_demux.csv",
G24="day12rep2_transcriptome/HTO_demux.csv",
G18="day14rep1_transcriptome/HTO_demux.csv",
G25="day14rep2_transcriptome/HTO_demux.csv",
G19="day19rep1_transcriptome/HTO_demux.csv",
G26="day19rep2_transcriptome/HTO_demux.csv",
G20="day28rep1_transcriptome/HTO_demux.csv",
G27="day28rep2_transcriptome/HTO_demux.csv"
)

demuxdatapathsnew=lapply(demuxfilenamesnew, function(x) paste0(countdataroots$new, x))

demuxdatapathsold=lapply(demuxfilenamesold, function(x) paste0(countdataroots$old, x))

demuxdatapaths= c(demuxdatapathsold, demuxdatapathsnew)

################################################################################
#map of days to datasets
################################################################################

dayrep.dataset.map=list(
  G1_GEX="G1",
  G2_GEX="G2",
  G3_GEX="G3",
  G4_GEX="G4",
  G5_GEX="G5",
  G6_GEX="G6",
  G7_GEX="G7",
  G8_GEX="G8",
  G9_GEX="G9",
  G10_GEX="G10",
  G11_GEX="G11",
  G12_GEX="G12",
  G13_GEX="G13",
day3rep1= "G14",
day3rep2= "G21",
day6rep1= "G15",
day6rep2= "G22",
day10rep1= "G16",
day10rep2= "G23",
day12rep1="G17",
day12rep2="G24",
day14rep1="G18",
day14rep2="G25",
day19rep1="G19",
day19rep2="G26",
day28rep1="G20",
day28rep2="G27"
)

################################################################################
#map replicate number directly
################################################################################
repmap=list(
G1= "R1",
G2= "R1",
G3= "R2",
G4= "R2",
G5= "R1",
G6= "R1",
G7= "R1",
G8= "R3",
G9= "R3",
G10="R3",
G11="R3",
G12="R3",
G13= "R2", 
G14= "R4",
G21= "R5",
G15= "R4",
G22= "R5",
G16= "R4",
G23= "R5",
G17="R4",
G24="R5",
G18="R4",
G25="R5",
G19="R4",
G26="R5",
G20="R4",
G27="R5"
)


fixoldcondition= Vectorize(function(x){
  gsub("x$", "1q", paste0("c", x))
  
}, USE.NAMES=F )

################################################################################
# fix the column names in a dataframe
################################################################################

fixnames= function(x){nms= x %>% colnames %>% make.names; colnames(x)=nms; x}

################################################################################
# sample barcode relationship tables
################################################################################

drfun=Vectorize(function(x){ (gsub("H7WT", "H7_WT", x) %>% strsplit(., split="_"))[[1]][c(3,4)] %>% paste(., collapse="")}, USE.NAMES=F)

sample.bc.rel.old=fread("~/metadata/samples_multiseq.csv") %>% as.data.frame %>% fixnames
sample.bc.rel.new=fread("~/metadata/BCsheet1_10x_annotationsheet_2023-04_ncnb.csv") %>% as.data.frame %>% fixnames %>% mutate(dayrep= drfun(Demultiplexed.Sample.Name), experiment_group= Vectorize(function(x) dayrep.dataset.map[[x]], USE.NAMES=F)(dayrep), sample_name=Demultiplexed.Sample.Name, multiseq_id=Multiplexing.Barcode.ID) 

barcode.relationships=list(
G1= sample.bc.rel.old %>% filter(experiment_group=="G1") %>%  select(experiment_group, sample_name, multiseq_id),
G2= sample.bc.rel.old %>% filter(experiment_group=="G2") %>%  select(experiment_group, sample_name, multiseq_id),
G3= sample.bc.rel.old %>% filter(experiment_group=="G3") %>%  select(experiment_group, sample_name, multiseq_id),
G4= sample.bc.rel.old %>% filter(experiment_group=="G4") %>%  select(experiment_group, sample_name, multiseq_id),
G5= sample.bc.rel.old %>% filter(experiment_group=="G5") %>%  select(experiment_group, sample_name, multiseq_id),
G6= sample.bc.rel.old %>% filter(experiment_group=="G6") %>%  select(experiment_group, sample_name, multiseq_id),
G7= sample.bc.rel.old %>% filter(experiment_group=="G7") %>%  select(experiment_group, sample_name, multiseq_id),
G8= sample.bc.rel.old %>% filter(experiment_group=="G8") %>%  select(experiment_group, sample_name, multiseq_id),
G9= sample.bc.rel.old %>% filter(experiment_group=="G9") %>%  select(experiment_group, sample_name, multiseq_id),
G10=sample.bc.rel.old %>% filter(experiment_group=="G10") %>% select(experiment_group, sample_name, multiseq_id),
G11=sample.bc.rel.old %>% filter(experiment_group=="G11") %>% select(experiment_group, sample_name, multiseq_id),
G12=sample.bc.rel.old %>% filter(experiment_group=="G12") %>% select(experiment_group, sample_name, multiseq_id),
G13=sample.bc.rel.old %>% filter(experiment_group=="G13") %>% select(experiment_group, sample_name, multiseq_id),
G14=sample.bc.rel.new %>% filter(experiment_group=="G14", grepl("CMO", sample_name)) %>% select(experiment_group, sample_name, multiseq_id),
G21=sample.bc.rel.new %>% filter(experiment_group=="G21", grepl("CMO", sample_name)) %>% select(experiment_group, sample_name, multiseq_id),
G15=sample.bc.rel.new %>% filter(experiment_group=="G15", grepl("CMO", sample_name)) %>% select(experiment_group, sample_name, multiseq_id),
G22=sample.bc.rel.new %>% filter(experiment_group=="G22", grepl("CMO", sample_name)) %>% select(experiment_group, sample_name, multiseq_id),
G16=sample.bc.rel.new %>% filter(experiment_group=="G16", grepl("CMO", sample_name)) %>% select(experiment_group, sample_name, multiseq_id),
G23=sample.bc.rel.new %>% filter(experiment_group=="G23", grepl("CMO", sample_name)) %>% select(experiment_group, sample_name, multiseq_id),
G17=sample.bc.rel.new %>% filter(experiment_group=="G17", grepl("CMO", sample_name)) %>% select(experiment_group, sample_name, multiseq_id),
G24=sample.bc.rel.new %>% filter(experiment_group=="G24", grepl("CMO", sample_name)) %>% select(experiment_group, sample_name, multiseq_id),
G18=sample.bc.rel.new %>% filter(experiment_group=="G18", grepl("CMO", sample_name)) %>% select(experiment_group, sample_name, multiseq_id),
G25=sample.bc.rel.new %>% filter(experiment_group=="G25", grepl("CMO", sample_name)) %>% select(experiment_group, sample_name, multiseq_id),
G19=sample.bc.rel.new %>% filter(experiment_group=="G19", grepl("CMO", sample_name)) %>% select(experiment_group, sample_name, multiseq_id),
G26=sample.bc.rel.new %>% filter(experiment_group=="G26", grepl("CMO", sample_name)) %>% select(experiment_group, sample_name, multiseq_id),
G20=sample.bc.rel.new %>% filter(experiment_group=="G20", grepl("CMO", sample_name)) %>% select(experiment_group, sample_name, multiseq_id),
G27=sample.bc.rel.new %>% filter(experiment_group=="G27", grepl("CMO", sample_name)) %>% select(experiment_group, sample_name, multiseq_id)
)  



  bcrels=barcode.relationships  %>% Reduce(rbind, .) %>% as.data.frame
  
  
  allds=bcrels %>% pull(experiment_group) %>% unique
  
  # incorporate unassigned, multiplet, blank categories
  
  bcrels=lapply(allds, function(x){
  
    rs= bcrels %>% filter(experiment_group==x)
    rbind(rs,
    rs[1,] %>% mutate (sample_name="Unassigned", multiseq_id="Unassigned"),
    rs[1, ] %>% mutate (sample_name="Multiplet", multiseq_id="Multiplet"),
    rs[1, ] %>% mutate (sample_name="Blank", multiseq_id="Blank")
    )
      
  }) %>% Reduce(rbind, .)
  
  
################################################################################
# neural crest cell type markers from Kameneva et al, 2021 
################################################################################
  
adamarkers=list(
  SCP=c("SOX10", "PLP1", "FOXD3", "FABP7", "S100B", "ERBB3", "NGFR", "MPZ", "COL2A1", "POSTN", "MOXD1", "GAS7", "ASCL1"),
  sympathoblasts=c("STMN2", "ELAVL4", "STMN4", "ISL1", "PRPH", "ELAVL2", "HMX1", "PHOX2B", "GATA3"),
  mesenchyme=c("COL1A1", "COL1A2", "COL12A1", "COL3A1", "VIM", "CXCL12", "TWIST1", "TWIST2"),
  intermediate_mesoderm=c( "GATA4", "HAND2"),
  melanocytes=c("MITF", "DCT"),
  endothelium=c("PECAM1", "KDR", "CAVIN2", "FLT1", "EFFL7", "PRCP"), 
  stemcells=c("POU5F1", "NANOG", "SOX2", "PODXL"),
  cortex=c("STAR", "NR5A1", "CYP17A1", "CYP11A1", "CYP21A2"), 
  erythroid= c("HBA2", "ALAS2"),
  kidney=c("PAX2", "LYPD1", "LHX1"),
  chromaffin=c("PNMT", "CHGA"),
  liver=c("HNF4A", "AHSG", "ITIH1")
)  
  

  
################################################################################
  # 17q and 1q centromere annotation ( for hg19)
################################################################################

  
centromeres= list(chr1=c(start=121484115, end=142582680), chr17=c(start=22262240, end=25270050)) %>% as.data.frame %>% t %>% as.data.frame %>% 
names2col(., coln="seqnames") %>% mutate(region="centromere") #%>% pivot_longer(., c("start", "end"), names_to="coordinate.type", values_to="position")

      
################################################################################
# make bam path
################################################################################

  
getbams=function(x) {

  truepath=paste0(config$out_root_host, "/scrna/",countdatafolders[[expbatches[[x]]]],"/",  dsname.dayrep.map[[x]], "/outs/per_sample_outs/")
  
  fls=list.files(paste0("~/mnt_out/scrna/", countdatafolders[[expbatches[[x]]]],"/", dsname.dayrep.map[[x]], "/outs/per_sample_outs/"))
  
  fls
  str1=paste0("mkdir  ~/path/to/bamfiles/",dsname.dayrep.map[[x]], "\n")
  
  
    str2=lapply(fls, function(fl) {
  
paste0(paste0(truepath, fl, "/count/sample_alignments.bam\n")) 
}) %>% Reduce(paste0,.)
cat(str2)
    
}
  


################################################################################
#import parse table data
################################################################################

gentable=fread(file.path(params$outpath, "parse_repliate_counts_per_day_mod.csv"))
colnames(gentable)=make.names(colnames(gentable))
gentable= gentable %>% mutate(actual.cgenotype=paste0("c", actual.genotype), stage=paste0("D", day))
gentable[["condition"]]=gsub( "cMYCN", "cWTMYCN", gentable[["condition"]])

genomap0=cbind(gentable[["condition"]], gentable[["actual.cgenotype"]])
genomap=genomap0[!duplicated(genomap0),]
get.actual.genotype=Vectorize(function(g){
  genomap[,2][genomap[,1]==g]
},USE.NAMES=T)






cellranger.config.files=list(
  
"/home/rstudio/mnt_out/metadata/cellrangerconfig_G1_GEX.csv",
"/home/rstudio/mnt_out/metadata/cellrangerconfig_G2_GEX.csv",
"/home/rstudio/mnt_out/metadata/cellrangerconfig_G3_GEX.csv",
"/home/rstudio/mnt_out/metadata/cellrangerconfig_G4_GEX.csv",
"/home/rstudio/mnt_out/metadata/cellrangerconfig_G5_GEX.csv",
"/home/rstudio/mnt_out/metadata/cellrangerconfig_G6_GEX.csv",
"/home/rstudio/mnt_out/metadata/cellrangerconfig_G7_GEX.csv",
"/home/rstudio/mnt_out/metadata/cellrangerconfig_G8_GEX.csv",
"/home/rstudio/mnt_out/metadata/cellrangerconfig_G9_GEX.csv",
"/home/rstudio/mnt_out/metadata/cellrangerconfig_G10_GEX.csv",
"/home/rstudio/mnt_out/metadata/cellrangerconfig_G11_GEX.csv",
"/home/rstudio/mnt_out/metadata/cellrangerconfig_G12_GEX.csv",
"/home/rstudio/mnt_out/metadata/cellrangerconfig_G13_GEX.csv",
"/home/rstudio/mnt_out/metadata/cellrangerconfig_G14_day3rep1.csv",
"/home/rstudio/mnt_out/metadata/cellrangerconfig_G15_day6rep1.csv",
"/home/rstudio/mnt_out/metadata/cellrangerconfig_G16_day10rep1.csv",
"/home/rstudio/mnt_out/metadata/cellrangerconfig_G17_day12rep1.csv",
"/home/rstudio/mnt_out/metadata/cellrangerconfig_G18_day14rep1.csv",
"/home/rstudio/mnt_out/metadata/cellrangerconfig_G19_day19rep1.csv",
"/home/rstudio/mnt_out/metadata/cellrangerconfig_G20_day28rep1.csv",
"/home/rstudio/mnt_out/metadata/cellrangerconfig_G21_day3rep2.csv",
"/home/rstudio/mnt_out/metadata/cellrangerconfig_G22_day6rep2.csv",
"/home/rstudio/mnt_out/metadata/cellrangerconfig_G23_day10rep2.csv",
"/home/rstudio/mnt_out/metadata/cellrangerconfig_G24_day12rep2.csv",
"/home/rstudio/mnt_out/metadata/cellrangerconfig_G25_day14rep2.csv",
"/home/rstudio/mnt_out/metadata/cellrangerconfig_G26_day19rep2.csv",
"/home/rstudio/mnt_out/metadata/cellrangerconfig_G27_day28rep2.csv"
)

names(cellranger.config.files)=paste0("G", 1:27)



    