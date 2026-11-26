/* CrisprGuard Interactive Web Laboratory Controller */

(() => {
  "use strict";

  const I18N = {
    en: {
      subtitle: "Cas9/Cas12a Guide Designer · CFD Off-Target Kinetics · Bae MMEJ Knockout Lab",
      exportCsv: "Export CSV",
      exportFasta: "Export FASTA",
      enzymeLabel: "CRISPR Endonuclease",
      presetLabel: "Curated Target Gene Presets",
      targetSeqLabel: "Target DNA / FASTA Sequence",
      bgSeqLabel: "Off-Target Background Reference (Multi-FASTA)",
      maxMmLabel: "Max Mismatches",
      nonCanonicalLabel: "Include non-canonical PAMs (NAG/NGA · TTTT/CTTV)",
      runAnalysis: "Run Design & Kinetics",
      customPreset: "Custom FASTA",
      customDesc: "Paste any DNA or FASTA target exon and an optional multi-FASTA background reference for CFD kinetics.",
      kpiTarget: "Target Length & GC",
      kpiGuides: "Candidate Guides",
      kpiTopScore: "Top Composite Score",
      kpiBestMmej: "Best MMEJ Knockout",
      kpiOffTargets: "Background Scanned",
      trackTitle: "Genomic Target Track & Guide RNA Map",
      trackHint: "Click any sgRNA/crRNA arrow to inspect exact DSB cleavage, MMEJ microhomology deletions & CFD off-target kinetics.",
      zoomLabel: "Zoom",
      tableTitle: "Ranked Candidate Guides",
      filterAll: "All",
      hidePol3: "Hide TTTT (Pol III)",
      searchPlaceholder: "Filter spacer / PAM...",
      colStrand: "Strand",
      colCut: "Cut bp",
      colSpacer: "Spacer (5'→3') + PAM",
      colOnTarget: "On-Target",
      colSpec: "Specificity",
      colMmej: "MMEJ OOF%",
      colComp: "Composite",
      colFlags: "Flags",
      copySpacer: "Copy Spacer",
      copied: "Copied!",
      spacerSeedLabel: "Protospacer Architecture & PAM-Proximal Seed",
      metricOnTarget: "On-Target Efficiency",
      metricSpec: "CFD Specificity",
      metricMmej: "MMEJ Out-of-Frame",
      onTargetBreakdownTitle: "On-Target Biophysical Decomposition",
      mmejDiagramTitle: "Microhomology (MMEJ) Deletion Alignment Around DSB Cut Site",
      offTargetMatrixTitle: "CFD Off-Target Mismatch Matrix (PAM-Proximal Seed Highlighted)",
      legSeed: "Seed MM",
      legDistal: "Distal MM",
      legWobble: "rG:dT Wobble",
      running: "Scanning PAM sites, CFD kinetics & MMEJ deletions...",
      doneStatus: (n, ms) => `${n} guides scored in ${ms} ms`,
      noGuides: "No candidate guides match the selected filters.",
      noOffTargets: "Zero off-target loci found within the mismatch threshold (100.0 specificity).",
      noMmej: "No 2–8 bp flanking microhomology repeats detected within ±25 bp of the DSB cut site.",
      frameshift: "Frameshift KO",
      inFrame: "In-Frame",
      cutBlunt: (pos) => `Blunt DSB @ bp ${pos}`,
      cutStaggered: (nts, ts) => `5-nt 5' Overhang (${nts} / ${ts})`,
      warnings: {
        pol3_terminator_tttt: "Pol III TTTT Terminator",
        strong_hairpin: "Strong Hairpin",
        moderate_hairpin: "Moderate Hairpin",
        scaffold_complementarity: "Scaffold Pairing",
        g_quadruplex_gggg: "G-Quadruplex (GGGG)",
        homopolymer_run: "Homopolymer >=4",
        extreme_gc: "Extreme GC%",
        seed_poly_u: "Seed Poly-U",
        high_risk_offtarget: "High-Risk Off-Target",
        in_frame_mmej_bias: "In-Frame MMEJ Bias",
      },
      compNames: {
        gc_balance: "GC Balance",
        position_prefs: "Position Prefs",
        seed_u_depletion: "Seed U-Depletion",
        pol3_terminator: "Pol III TTTT",
        homopolymer: "Homopolymer",
        secondary_structure: "Hairpin / Fold",
        tm_stability: "Tm Stability",
      },
    },
    tr: {
      subtitle: "Cas9/Cas12a Kılavuz Tasarımı · CFD Hedef-Dışı Kinetiği · Bae MMEJ Nakavt Labı",
      exportCsv: "CSV İndir",
      exportFasta: "FASTA İndir",
      enzymeLabel: "CRISPR Endonükleaz",
      presetLabel: "Küratörlü Hedef Gen Örnekleri",
      targetSeqLabel: "Hedef DNA / FASTA Dizisi",
      bgSeqLabel: "Hedef-Dışı Arka Plan Referansı (Multi-FASTA)",
      maxMmLabel: "Maks. Uyumsuzluk",
      nonCanonicalLabel: "Kanonik dışı PAM'ları dahil et (NAG/NGA · TTTT/CTTV)",
      runAnalysis: "Tasarım ve Kinetiği Çalıştır",
      customPreset: "Özel FASTA",
      customDesc: "Kendi DNA veya FASTA hedef ekzonunuzu ve isteğe bağlı arka plan referans dizilerinizi yapıştırın.",
      kpiTarget: "Hedef Uzunluk & GC",
      kpiGuides: "Aday Kılavuzlar",
      kpiTopScore: "En Yüksek Bileşik Skor",
      kpiBestMmej: "En İyi MMEJ Nakavt",
      kpiOffTargets: "Taranan Arka Plan",
      trackTitle: "Genomik Hedef İzi ve Kılavuz RNA Haritası",
      trackHint: "Kesin DSB kesim bölgesini, MMEJ mikrohomoloji delesyonlarını ve CFD hedef-dışı matrisini incelemek için herhangi bir oka tıklayın.",
      zoomLabel: "Yakınlaştır",
      tableTitle: "Sıralanmış Aday Kılavuz RNA'lar",
      filterAll: "Tümü",
      hidePol3: "TTTT (Pol III) Gizle",
      searchPlaceholder: "Spacer / PAM filtrele...",
      colStrand: "İplik",
      colCut: "Kesim bp",
      colSpacer: "Spacer (5'→3') + PAM",
      colOnTarget: "Hedef-İçi",
      colSpec: "Özgüllük",
      colMmej: "MMEJ OOF%",
      colComp: "Bileşik",
      colFlags: "Uyarılar",
      copySpacer: "Spacer Kopyala",
      copied: "Kopyalandı!",
      spacerSeedLabel: "Protospacer Mimarisi ve PAM-Yakını Çekirdek (Seed)",
      metricOnTarget: "Hedef-İçi Etkinlik",
      metricSpec: "CFD Özgüllüğü",
      metricMmej: "MMEJ Çerçeve Kayması",
      onTargetBreakdownTitle: "Hedef-İçi Biyofiziksel Bileşenler",
      mmejDiagramTitle: "DSB Kesim Noktası Çevresi Mikrohomoloji (MMEJ) Delesyon Hizalaması",
      offTargetMatrixTitle: "CFD Hedef-Dışı Uyumsuzluk Matrisi (PAM-Yakını Çekirdek Vurgulu)",
      legSeed: "Çekirdek MM",
      legDistal: "Uzak MM",
      legWobble: "rG:dT Wobble",
      running: "PAM bölgeleri, CFD kinetiği ve MMEJ delesyonları taranıyor...",
      doneStatus: (n, ms) => `${n} kılavuz ${ms} ms içinde puanlandı`,
      noGuides: "Seçilen filtrelere uyan aday kılavuz bulunamadı.",
      noOffTargets: "Belirtilen uyumsuzluk eşiği içinde hedef-dışı lokus bulunamadı (100.0 özgüllük).",
      noMmej: "DSB kesim noktasının ±25 bp çevresinde 2–8 bp mikrohomoloji tekrarı saptanmadı.",
      frameshift: "Çerçeve Kayması (KO)",
      inFrame: "Çerçeve İçi",
      cutBlunt: (pos) => `Küt DSB @ bp ${pos}`,
      cutStaggered: (nts, ts) => `5-nt 5' Çıkıntı (${nts} / ${ts})`,
      warnings: {
        pol3_terminator_tttt: "Pol III TTTT Sonlandırıcı",
        strong_hairpin: "Güçlü Saç Tokası",
        moderate_hairpin: "Orta Saç Tokası",
        scaffold_complementarity: "İskele Eşleşmesi",
        g_quadruplex_gggg: "G-Kuadrupleks (GGGG)",
        homopolymer_run: "Homopolimer >=4",
        extreme_gc: "Uç GC%",
        seed_poly_u: "Çekirdek Poli-U",
        high_risk_offtarget: "Yüksek Riskli Hedef-Dışı",
        in_frame_mmej_bias: "Çerçeve-İçi MMEJ Eğilimi",
      },
      compNames: {
        gc_balance: "GC Dengesi",
        position_prefs: "Pozisyon Tercihi",
        seed_u_depletion: "Çekirdek U-Azlığı",
        pol3_terminator: "Pol III TTTT",
        homopolymer: "Homopolimer",
        secondary_structure: "İkincil Yapı",
        tm_stability: "Tm Kararlılığı",
      },
    },
  };

  const state = {
    lang: "en",
    enzyme: "spcas9",
    samples: [],
    activeSampleId: "tp53_exon5",
    report: null,
    filteredGuides: [],
    selectedGuideId: null,
    hoveredGuideId: null,
    strandFilter: "all",
    hidePol3: false,
    searchQuery: "",
    zoom: 1.0,
    hitBoxes: [],
    activeMmeJIdx: 0,
  };

  const els = {
    langEnBtn: document.getElementById("langEnBtn"),
    langTrBtn: document.getElementById("langTrBtn"),
    exportCsvBtn: document.getElementById("exportCsvBtn"),
    exportFastaBtn: document.getElementById("exportFastaBtn"),
    enzymeToggle: document.getElementById("enzymeToggle"),
    presetPills: document.getElementById("presetPills"),
    presetDesc: document.getElementById("presetDesc"),
    targetInput: document.getElementById("targetInput"),
    bgInput: document.getElementById("bgInput"),
    targetCharCount: document.getElementById("targetCharCount"),
    bgCharCount: document.getElementById("bgCharCount"),
    maxMmSelect: document.getElementById("maxMmSelect"),
    nonCanonicalPamCheck: document.getElementById("nonCanonicalPamCheck"),
    runBtn: document.getElementById("runBtn"),
    statusMessage: document.getElementById("statusMessage"),
    kpiTargetVal: document.getElementById("kpiTargetVal"),
    kpiTargetSub: document.getElementById("kpiTargetSub"),
    kpiGuidesVal: document.getElementById("kpiGuidesVal"),
    kpiGuidesSub: document.getElementById("kpiGuidesSub"),
    kpiTopVal: document.getElementById("kpiTopVal"),
    kpiTopSub: document.getElementById("kpiTopSub"),
    kpiMmejVal: document.getElementById("kpiMmejVal"),
    kpiMmejSub: document.getElementById("kpiMmejSub"),
    kpiBgVal: document.getElementById("kpiBgVal"),
    kpiBgSub: document.getElementById("kpiBgSub"),
    canvasViewport: document.getElementById("canvasViewport"),
    canvas: document.getElementById("genomicCanvas"),
    canvasTooltip: document.getElementById("canvasTooltip"),
    zoomSlider: document.getElementById("zoomSlider"),
    zoomVal: document.getElementById("zoomVal"),
    tableCountSub: document.getElementById("tableCountSub"),
    strandFilter: document.getElementById("strandFilter"),
    hidePol3Check: document.getElementById("hidePol3Check"),
    guideSearchInput: document.getElementById("guideSearchInput"),
    guideTableBody: document.getElementById("guideTableBody"),
    drawerBackdrop: document.getElementById("drawerBackdrop"),
    inspectorDrawer: document.getElementById("inspectorDrawer"),
    closeDrawerBtn: document.getElementById("closeDrawerBtn"),
    copySpacerBtn: document.getElementById("copySpacerBtn"),
    drawerGuideId: document.getElementById("drawerGuideId"),
    drawerTitle: document.getElementById("drawerTitle"),
    drawerRankBadge: document.getElementById("drawerRankBadge"),
    drawerMeta: document.getElementById("drawerMeta"),
    drawerCutInfo: document.getElementById("drawerCutInfo"),
    drawerSpacerViz: document.getElementById("drawerSpacerViz"),
    drawerWarnings: document.getElementById("drawerWarnings"),
    gaugeOnTarget: document.getElementById("gaugeOnTarget"),
    gaugeSpec: document.getElementById("gaugeSpec"),
    gaugeMmej: document.getElementById("gaugeMmej"),
    valOnTarget: document.getElementById("valOnTarget"),
    valSpec: document.getElementById("valSpec"),
    valMmej: document.getElementById("valMmej"),
    gradeMmej: document.getElementById("gradeMmej"),
    subOnTarget: document.getElementById("subOnTarget"),
    subSpec: document.getElementById("subSpec"),
    subMmej: document.getElementById("subMmej"),
    onTargetMetaGrid: document.getElementById("onTargetMetaGrid"),
    onTargetComponents: document.getElementById("onTargetComponents"),
    mmejRatioText: document.getElementById("mmejRatioText"),
    mmejFsBar: document.getElementById("mmejFsBar"),
    mmejIfBar: document.getElementById("mmejIfBar"),
    mmejRefTrack: document.getElementById("mmejRefTrack"),
    mmejPatternsList: document.getElementById("mmejPatternsList"),
    offTargetMatrix: document.getElementById("offTargetMatrix"),
  };

  const t = (key) => I18N[state.lang][key] ?? I18N.en[key] ?? key;

  function escapeHtml(str) {
    return String(str)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;");
  }

  function countBases(text) {
    return text
      .split(/\r?\n/)
      .filter((line) => !line.trim().startsWith(">") && !line.trim().startsWith(";"))
      .join("")
      .replace(/[^ACGTUNacgtun]/g, "").length;
  }

  function updateCharCounts() {
    els.targetCharCount.textContent = `${countBases(els.targetInput.value)} bp`;
    els.bgCharCount.textContent = `${countBases(els.bgInput.value)} bp`;
  }

  function scoreColor(score) {
    if (score >= 75) return "#10b981";
    if (score >= 60) return "#06b6d4";
    if (score >= 45) return "#f59e0b";
    return "#f43f5e";
  }

  function scorePillClass(score) {
    if (score >= 75) return "pill-emerald";
    if (score >= 60) return "pill-cyan";
    if (score >= 45) return "pill-amber";
    return "pill-rose";
  }

  function applyLanguage() {
    document.documentElement.lang = state.lang;
    els.langEnBtn.classList.toggle("active", state.lang === "en");
    els.langTrBtn.classList.toggle("active", state.lang === "tr");

    document.querySelectorAll("[data-i18n]").forEach((node) => {
      const key = node.getAttribute("data-i18n");
      const val = t(key);
      if (typeof val === "string") node.textContent = val;
    });
    document.querySelectorAll("[data-i18n-placeholder]").forEach((node) => {
      const key = node.getAttribute("data-i18n-placeholder");
      const val = t(key);
      if (typeof val === "string") node.setAttribute("placeholder", val);
    });

    renderPresets();
    updatePresetDescription();
    if (state.report) {
      renderTable();
      drawCanvas();
      if (!els.inspectorDrawer.classList.contains("hidden") && state.selectedGuideId) {
        const g = state.report.guides.find((x) => x.id === state.selectedGuideId);
        if (g) renderInspector(g);
      }
    }
  }

  function renderPresets() {
    els.presetPills.innerHTML = "";
    state.samples.forEach((sample) => {
      const btn = document.createElement("button");
      btn.type = "button";
      btn.className = `preset-pill ${state.activeSampleId === sample.id ? "active" : ""}`;
      btn.textContent = sample.name;
      btn.addEventListener("click", () => selectSample(sample.id, true));
      els.presetPills.appendChild(btn);
    });
    const customBtn = document.createElement("button");
    customBtn.type = "button";
    customBtn.className = `preset-pill ${state.activeSampleId === "custom" ? "active" : ""}`;
    customBtn.textContent = t("customPreset");
    customBtn.addEventListener("click", () => {
      state.activeSampleId = "custom";
      renderPresets();
      updatePresetDescription();
      els.targetInput.focus();
    });
    els.presetPills.appendChild(customBtn);
  }

  function updatePresetDescription() {
    if (state.activeSampleId === "custom") {
      els.presetDesc.textContent = t("customDesc");
      return;
    }
    const sample = state.samples.find((s) => s.id === state.activeSampleId);
    if (!sample) return;
    const desc = state.lang === "tr" ? sample.description_tr : sample.description_en;
    els.presetDesc.textContent = `${sample.locus} — ${desc}`;
  }

  function selectSample(sampleId, autoRun = false) {
    const sample = state.samples.find((s) => s.id === sampleId);
    if (!sample) return;
    state.activeSampleId = sample.id;
    els.targetInput.value = sample.target_fasta.trim();
    els.bgInput.value = sample.background_fasta.trim();
    updateCharCounts();
    renderPresets();
    updatePresetDescription();
    if (autoRun) runDesign();
  }

  async function runDesign() {
    const seq = els.targetInput.value.trim();
    if (!seq) {
      els.statusMessage.textContent = "Target sequence cannot be empty.";
      els.statusMessage.classList.add("error");
      return;
    }
    els.runBtn.disabled = true;
    els.statusMessage.classList.remove("error");
    els.statusMessage.textContent = t("running");
    const t0 = performance.now();

    try {
      const resp = await fetch("/api/design", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          sequence: seq,
          background: els.bgInput.value.trim() || null,
          enzyme: state.enzyme,
          max_mismatches: parseInt(els.maxMmSelect.value, 10),
          include_noncanonical_pam: els.nonCanonicalPamCheck.checked,
        }),
      });
      const data = await resp.json();
      if (!resp.ok) {
        throw new Error(data.error || `HTTP ${resp.status}`);
      }
      state.report = data;
      state.selectedGuideId = data.guides.length > 0 ? data.guides[0].id : null;
      state.activeMmeJIdx = 0;

      const elapsed = Math.max(1, Math.round(performance.now() - t0));
      els.statusMessage.textContent = t("doneStatus")(data.total_candidates, elapsed);
      els.exportCsvBtn.disabled = data.guides.length === 0;
      els.exportFastaBtn.disabled = data.guides.length === 0;

      renderKPIs(data);
      applyFiltersAndRender();
      drawCanvas();
    } catch (err) {
      els.statusMessage.textContent = err.message;
      els.statusMessage.classList.add("error");
    } finally {
      els.runBtn.disabled = false;
    }
  }

  function renderKPIs(report) {
    els.kpiTargetVal.textContent = `${report.target_length} bp`;
    els.kpiTargetSub.textContent = `GC ${report.target_gc_percent}% · ${report.target_name}`;

    els.kpiGuidesVal.textContent = `${report.total_candidates}`;
    els.kpiGuidesSub.textContent = `+${report.plus_strand_count} / -${report.minus_strand_count} (${report.enzyme_spec.pam_motif})`;

    if (report.guides.length > 0) {
      const best = report.guides[0];
      els.kpiTopVal.textContent = `${best.composite_score.toFixed(1)}`;
      els.kpiTopSub.textContent = `${best.id} (${best.strand}) · OnTgt ${best.on_target_score} · Spec ${best.specificity_score}`;

      const bestKo = [...report.guides].sort(
        (a, b) => b.mmej.out_of_frame_probability - a.mmej.out_of_frame_probability
      )[0];
      els.kpiMmejVal.textContent = `${bestKo.mmej.out_of_frame_probability.toFixed(1)}%`;
      els.kpiMmejSub.textContent = `${bestKo.id} · ${bestKo.mmej.knockout_grade}`;
    } else {
      els.kpiTopVal.textContent = "—";
      els.kpiTopSub.textContent = "0 guides";
      els.kpiMmejVal.textContent = "—";
      els.kpiMmejSub.textContent = "—";
    }

    els.kpiBgVal.textContent = `${report.background_bp.toLocaleString()} bp`;
    els.kpiBgSub.textContent = `${report.background_contigs} contigs · <=${report.max_mismatches} mm`;
  }

  function applyFiltersAndRender() {
    if (!state.report) return;
    const q = state.searchQuery.trim().toUpperCase();
    state.filteredGuides = state.report.guides.filter((g) => {
      if (state.strandFilter !== "all" && g.strand !== state.strandFilter) return false;
      if (state.hidePol3 && g.on_target_detail.pol3_terminator) return false;
      if (q) {
        const hay = `${g.id} ${g.spacer} ${g.pam}`.toUpperCase();
        if (!hay.includes(q)) return false;
      }
      return true;
    });
    renderTable();
    drawCanvas();
  }

  function formatSpacerWithSeed(guide) {
    const sp = guide.spacer;
    if (guide.enzyme === "spcas9") {
      // Positions 1..10 distal, 11..20 seed
      return `${escapeHtml(sp.slice(0, 10))}<span class="seed-underline">${escapeHtml(sp.slice(10))}</span><span class="pam-chip">${escapeHtml(guide.pam)}</span>`;
    }
    // LbCas12a: 5' PAM + positions 1..8 seed + 9..23 distal
    return `<span class="pam-chip" style="margin-left:0;margin-right:5px;">${escapeHtml(guide.pam)}</span><span class="seed-underline">${escapeHtml(sp.slice(0, 8))}</span>${escapeHtml(sp.slice(8))}`;
  }

  function renderTable() {
    const guides = state.filteredGuides;
    const total = state.report ? state.report.total_candidates : 0;
    els.tableCountSub.textContent = `${guides.length} / ${total}`;

    if (guides.length === 0) {
      els.guideTableBody.innerHTML = `<tr><td colspan="11" class="empty-note">${escapeHtml(t("noGuides"))}</td></tr>`;
      return;
    }

    const warnLabels = I18N[state.lang].warnings;
    els.guideTableBody.innerHTML = guides
      .map((g) => {
        const sel = g.id === state.selectedGuideId ? "selected" : "";
        const flagsHtml =
          g.warnings.length === 0
            ? `<span class="flag-ok">OK</span>`
            : g.warnings
                .slice(0, 2)
                .map((w) => `<span class="flag-chip">${escapeHtml(warnLabels[w] || w)}</span>`)
                .join("") + (g.warnings.length > 2 ? `<span class="mono-muted">+${g.warnings.length - 2}</span>` : "");

        return `
          <tr class="${sel}" data-guide-id="${escapeHtml(g.id)}">
            <td class="mono-muted">${g.rank}</td>
            <td><span class="id-badge">${escapeHtml(g.id)}</span></td>
            <td><span class="strand-tag ${g.strand === "+" ? "strand-plus" : "strand-minus"}">${g.strand}</span></td>
            <td class="mono-muted">${g.cut_pos}</td>
            <td class="seq-cell">${formatSpacerWithSeed(g)}</td>
            <td class="mono-muted">${g.gc_percent.toFixed(1)}%</td>
            <td>
              <div class="bar-cell">
                <div class="mini-bar"><span style="width:${g.on_target_score}%;background:${scoreColor(g.on_target_score)}"></span></div>
                <span>${g.on_target_score.toFixed(1)}</span>
              </div>
            </td>
            <td>
              <div class="bar-cell">
                <div class="mini-bar"><span style="width:${g.specificity_score}%;background:#8b5cf6"></span></div>
                <span>${g.specificity_score.toFixed(1)}</span>
              </div>
            </td>
            <td>
              <div class="bar-cell">
                <div class="mini-bar"><span style="width:${g.mmej.out_of_frame_probability}%;background:#06b6d4"></span></div>
                <span>${g.mmej.out_of_frame_probability.toFixed(1)}%</span>
              </div>
            </td>
            <td><span class="score-pill ${scorePillClass(g.composite_score)}">${g.composite_score.toFixed(1)}</span></td>
            <td>${flagsHtml}</td>
          </tr>
        `;
      })
      .join("");

    els.guideTableBody.querySelectorAll("tr[data-guide-id]").forEach((row) => {
      row.addEventListener("click", () => {
        const gid = row.getAttribute("data-guide-id");
        selectGuide(gid, true);
      });
    });
  }

  /* ---------------- Canvas 2D Genomic Target Track ---------------- */
  function packLanes(guides) {
    const sorted = [...guides].sort((a, b) => a.start - b.start || a.end - b.end);
    const laneEnds = [];
    const placed = [];
    for (const g of sorted) {
      let laneIdx = -1;
      for (let i = 0; i < laneEnds.length; i++) {
        if (g.start >= laneEnds[i] + 3) {
          laneIdx = i;
          break;
        }
      }
      if (laneIdx === -1) {
        laneIdx = laneEnds.length;
        laneEnds.push(g.end);
      } else {
        laneEnds[laneIdx] = g.end;
      }
      placed.push({ guide: g, lane: laneIdx });
    }
    return { placed, laneCount: Math.max(1, laneEnds.length) };
  }

  function complementBase(b) {
    return { A: "T", T: "A", C: "G", G: "C" }[b] || "N";
  }

  function baseColor(b) {
    return { A: "#34d399", T: "#fb7185", C: "#38bdf8", G: "#fbbf24" }[b] || "#94a3b8";
  }

  function drawCanvas() {
    if (!state.report) return;
    const report = state.report;
    const guides = state.filteredGuides;
    const plusPack = packLanes(guides.filter((g) => g.strand === "+"));
    const minusPack = packLanes(guides.filter((g) => g.strand === "-"));

    const padX = 48;
    const laneH = 22;
    const laneGap = 5;
    const topLanesH = Math.max(2, Math.min(plusPack.laneCount, 8)) * (laneH + laneGap) + 28;
    const botLanesH = Math.max(2, Math.min(minusPack.laneCount, 8)) * (laneH + laneGap) + 24;
    const duplexH = 54;
    const rulerH = 30;
    const cssHeight = rulerH + topLanesH + duplexH + botLanesH;

    const viewportWidth = Math.max(760, els.canvasViewport.clientWidth - 2);
    const cssWidth = Math.round(viewportWidth * state.zoom);
    const dpr = window.devicePixelRatio || 1;

    els.canvas.style.width = `${cssWidth}px`;
    els.canvas.style.height = `${cssHeight}px`;
    els.canvas.width = Math.round(cssWidth * dpr);
    els.canvas.height = Math.round(cssHeight * dpr);

    const ctx = els.canvas.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

    // Clear canvas
    ctx.fillStyle = "#070c17";
    ctx.fillRect(0, 0, cssWidth, cssHeight);

    const seqLen = Math.max(1, report.target_length);
    const trackW = cssWidth - padX * 2;
    const bpToX = (bp) => padX + (bp / seqLen) * trackW;
    const bpWidth = trackW / seqLen;

    const duplexCenterY = rulerH + topLanesH + duplexH / 2;
    const topStrandY = duplexCenterY - 12;
    const botStrandY = duplexCenterY + 12;

    // 1. Coordinate Ruler & Grid
    const tickStep = seqLen <= 120 ? 10 : seqLen <= 350 ? 25 : 50;
    ctx.font = "600 10px 'JetBrains Mono', monospace";
    ctx.textAlign = "center";
    for (let bp = 0; bp <= seqLen; bp += tickStep) {
      const x = bpToX(bp);
      ctx.strokeStyle = "rgba(30, 48, 82, 0.45)";
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(x, rulerH);
      ctx.lineTo(x, cssHeight - 10);
      ctx.stroke();

      ctx.fillStyle = "#6b82a8";
      ctx.fillText(`${bp} bp`, x, 18);
    }

    // 2. Double-Stranded DNA Duplex
    ctx.strokeStyle = "#1e3a5f";
    ctx.lineWidth = 3;
    ctx.beginPath();
    ctx.moveTo(padX, topStrandY);
    ctx.lineTo(padX + trackW, topStrandY);
    ctx.moveTo(padX, botStrandY);
    ctx.lineTo(padX + trackW, botStrandY);
    ctx.stroke();

    // 5' / 3' labels
    ctx.font = "700 10px 'JetBrains Mono', monospace";
    ctx.fillStyle = "#34d399";
    ctx.textAlign = "right";
    ctx.fillText("5'", padX - 10, topStrandY + 3);
    ctx.fillStyle = "#a78bfa";
    ctx.fillText("3'", padX - 10, botStrandY + 3);
    ctx.textAlign = "left";
    ctx.fillStyle = "#34d399";
    ctx.fillText("3'", padX + trackW + 10, topStrandY + 3);
    ctx.fillStyle = "#a78bfa";
    ctx.fillText("5'", padX + trackW + 10, botStrandY + 3);

    // Base-pair rungs & nucleotide letters when zoomed in
    const rungStep = bpWidth >= 7 ? 1 : bpWidth >= 3.5 ? 2 : 4;
    ctx.lineWidth = 1;
    for (let i = 0; i < seqLen; i += rungStep) {
      const x = bpToX(i + 0.5);
      ctx.strokeStyle = "rgba(56, 189, 248, 0.22)";
      ctx.beginPath();
      ctx.moveTo(x, topStrandY + 2);
      ctx.lineTo(x, botStrandY - 2);
      ctx.stroke();

      if (bpWidth >= 8.5) {
        const b = report.target_sequence[i];
        ctx.font = "700 9px 'JetBrains Mono', monospace";
        ctx.textAlign = "center";
        ctx.fillStyle = baseColor(b);
        ctx.fillText(b, x, topStrandY - 4);
        ctx.fillStyle = baseColor(complementBase(b));
        ctx.fillText(complementBase(b), x, botStrandY + 11);
      }
    }

    // 3. Draw Guide Arrows on + and - strands
    state.hitBoxes = [];

    function drawGuideItem(item, isPlus) {
      const g = item.guide;
      const lane = item.lane % 8;
      const y = isPlus
        ? topStrandY - 24 - lane * (laneH + laneGap) - laneH
        : botStrandY + 22 + lane * (laneH + laneGap);

      const x1 = bpToX(g.start);
      const x2 = bpToX(g.end);
      const w = Math.max(12, x2 - x1);
      const h = 16;
      const col = scoreColor(g.composite_score);
      const isSelected = g.id === state.selectedGuideId;
      const isHovered = g.id === state.hoveredGuideId;

      state.hitBoxes.push({ guide: g, x: x1, y, w, h });

      // Spacer body
      ctx.save();
      if (isSelected || isHovered) {
        ctx.shadowColor = isSelected ? "#10b981" : "#06b6d4";
        ctx.shadowBlur = isSelected ? 12 : 8;
      }

      // Rounded bar
      ctx.fillStyle = isSelected ? col : `${col}cc`;
      ctx.beginPath();
      ctx.roundRect(x1, y, w, h, 4);
      ctx.fill();

      // PAM cap highlight
      const px1 = bpToX(g.pam_start);
      const px2 = bpToX(g.pam_end);
      const pw = Math.max(5, px2 - px1);
      ctx.fillStyle = "#8b5cf6";
      ctx.beginPath();
      ctx.roundRect(px1, y, pw, h, 3);
      ctx.fill();

      // Selection / Hover border
      if (isSelected || isHovered) {
        ctx.strokeStyle = "#f8fafc";
        ctx.lineWidth = isSelected ? 2 : 1.2;
        ctx.beginPath();
        ctx.roundRect(x1 - 0.5, y - 0.5, w + 1, h + 1, 4);
        ctx.stroke();
      }
      ctx.restore();

      // Pol III warning tick
      if (g.on_target_detail.pol3_terminator) {
        ctx.fillStyle = "#f43f5e";
        ctx.beginPath();
        ctx.arc(x1 + 5, y + h / 2, 3, 0, Math.PI * 2);
        ctx.fill();
      }

      // Label inside bar if wide enough
      if (w >= 46) {
        ctx.font = "700 9.5px 'JetBrains Mono', monospace";
        ctx.fillStyle = "#04111d";
        ctx.textAlign = "center";
        const label = `${g.id} ${isPlus ? "→" : "←"}`;
        ctx.fillText(label, x1 + w / 2, y + 11.5);
      }
    }

    plusPack.placed.forEach((item) => drawGuideItem(item, true));
    minusPack.placed.forEach((item) => drawGuideItem(item, false));

    // 4. Draw exact DSB Cut Site Indicator for the selected guide
    const selGuide = report.guides.find((g) => g.id === state.selectedGuideId);
    if (selGuide) {
      const cutX = bpToX(selGuide.cut_pos);
      ctx.save();
      ctx.strokeStyle = "#f43f5e";
      ctx.lineWidth = 2;
      ctx.setLineDash([4, 3]);
      ctx.beginPath();
      ctx.moveTo(cutX, rulerH + 6);
      ctx.lineTo(cutX, cssHeight - 10);
      ctx.stroke();

      // Staggered Cas12a target-strand cut line
      if (selGuide.cut_ts !== selGuide.cut_pos) {
        const tsX = bpToX(selGuide.cut_ts);
        ctx.strokeStyle = "#fb923c";
        ctx.beginPath();
        ctx.moveTo(tsX, duplexCenterY - 26);
        ctx.lineTo(tsX, duplexCenterY + 26);
        ctx.stroke();
      }

      ctx.setLineDash([]);
      ctx.fillStyle = "#f43f5e";
      ctx.font = "700 10px 'JetBrains Mono', monospace";
      ctx.textAlign = "center";
      ctx.fillText(`✂ ${selGuide.id} @ ${selGuide.cut_pos}`, cutX, rulerH + 14);
      ctx.restore();
    }
  }

  function findHitAt(clientX, clientY) {
    const rect = els.canvas.getBoundingClientRect();
    const x = clientX - rect.left;
    const y = clientY - rect.top;
    for (let i = state.hitBoxes.length - 1; i >= 0; i--) {
      const box = state.hitBoxes[i];
      if (x >= box.x - 2 && x <= box.x + box.w + 2 && y >= box.y - 2 && y <= box.y + box.h + 2) {
        return box.guide;
      }
    }
    return null;
  }

  function handleCanvasMove(ev) {
    const hit = findHitAt(ev.clientX, ev.clientY);
    const prev = state.hoveredGuideId;
    state.hoveredGuideId = hit ? hit.id : null;
    if (prev !== state.hoveredGuideId) drawCanvas();

    if (!hit) {
      els.canvasTooltip.classList.add("hidden");
      return;
    }

    els.canvasTooltip.innerHTML = `
      <div class="tt-head">
        <span>${escapeHtml(hit.id)} (#${hit.rank} · ${hit.strand})</span>
        <span style="color:${scoreColor(hit.composite_score)}">${hit.composite_score.toFixed(1)}</span>
      </div>
      <div class="tt-seq">${escapeHtml(hit.spacer)} <span style="color:#c4b5fd">${escapeHtml(hit.pam)}</span></div>
      <div class="tt-grid">
        <span>On-Target: <b>${hit.on_target_score.toFixed(1)}</b></span>
        <span>CFD Spec: <b>${hit.specificity_score.toFixed(1)}</b></span>
        <span>MMEJ OOF: <b>${hit.mmej.out_of_frame_probability.toFixed(1)}%</b></span>
        <span>Cut bp: <b>${hit.cut_pos}</b></span>
      </div>
    `;
    const left = Math.min(window.innerWidth - 260, ev.clientX + 14);
    const top = Math.min(window.innerHeight - 130, ev.clientY + 14);
    els.canvasTooltip.style.left = `${left}px`;
    els.canvasTooltip.style.top = `${top}px`;
    els.canvasTooltip.classList.remove("hidden");
  }

  function handleCanvasClick(ev) {
    const hit = findHitAt(ev.clientX, ev.clientY);
    if (hit) {
      selectGuide(hit.id, true);
    }
  }

  /* ---------------- Guide Inspector Drawer ---------------- */
  function selectGuide(guideId, openDrawer = true) {
    state.selectedGuideId = guideId;
    state.activeMmeJIdx = 0;
    renderTable();
    drawCanvas();
    if (!state.report) return;
    const guide = state.report.guides.find((g) => g.id === guideId);
    if (guide && openDrawer) {
      renderInspector(guide);
      els.drawerBackdrop.classList.remove("hidden");
      els.inspectorDrawer.classList.remove("hidden");
    }
  }

  function closeDrawer() {
    els.drawerBackdrop.classList.add("hidden");
    els.inspectorDrawer.classList.add("hidden");
  }

  function setGauge(circleEl, percent) {
    const circumference = 2 * Math.PI * 36;
    const clamped = Math.max(0, Math.min(100, percent));
    circleEl.style.strokeDasharray = `${circumference}`;
    circleEl.style.strokeDashoffset = `${circumference * (1 - clamped / 100)}`;
  }

  function renderInspector(g) {
    const isCas9 = g.enzyme === "spcas9";
    els.drawerGuideId.textContent = g.id;
    els.drawerRankBadge.textContent = `#${g.rank} · ${g.composite_score.toFixed(1)}`;
    els.drawerTitle.textContent = isCas9 ? "SpCas9 sgRNA" : "LbCas12a crRNA";
    els.drawerMeta.textContent = `${g.strand} strand · bp ${g.start + 1}–${g.end} · GC ${g.gc_percent}%`;
    els.drawerCutInfo.textContent = isCas9
      ? t("cutBlunt")(g.cut_pos)
      : t("cutStaggered")(g.cut_pos, g.cut_ts);

    // 1. Protospacer & PAM tile architecture with DSB cut bar
    const tiles = [];
    if (!isCas9) {
      for (let i = 0; i < g.pam.length; i++) {
        tiles.push(`<span class="nt-tile is-pam">${g.pam[i]}<small>P${i + 1}</small></span>`);
      }
    }
    for (let i = 0; i < g.spacer.length; i++) {
      const pos = i + 1;
      const inSeed = isCas9 ? pos >= 11 : pos <= 8;
      if (isCas9 && pos === 18) {
        tiles.push(`<span class="nt-cut-marker" title="DSB Cut (-3 bp)"></span>`);
      }
      if (!isCas9 && pos === 19) {
        tiles.push(`<span class="nt-cut-marker" title="NTS Cut (+18 nt)"></span>`);
      }
      tiles.push(
        `<span class="nt-tile ${inSeed ? "is-seed" : ""}">${g.spacer[i]}<small>${pos}</small></span>`
      );
    }
    if (isCas9) {
      for (let i = 0; i < g.pam.length; i++) {
        tiles.push(`<span class="nt-tile is-pam">${g.pam[i]}<small>P${i + 1}</small></span>`);
      }
    }
    els.drawerSpacerViz.innerHTML = tiles.join("");

    const warnMap = I18N[state.lang].warnings;
    els.drawerWarnings.innerHTML =
      g.warnings.length === 0
        ? `<span class="flag-ok">✓ No biophysical or Pol III risk flags</span>`
        : g.warnings
            .map((w) => `<span class="flag-chip">${escapeHtml(warnMap[w] || w)}</span>`)
            .join("");

    // 2. Metric Gauges
    els.valOnTarget.textContent = g.on_target_score.toFixed(1);
    els.subOnTarget.textContent = `Tm ${g.on_target_detail.tm_celsius}°C · GC ${g.gc_percent}%`;
    setGauge(els.gaugeOnTarget, g.on_target_score);

    const ot = g.offtarget_summary;
    els.valSpec.textContent = g.specificity_score.toFixed(1);
    els.subSpec.textContent = ot
      ? `${ot.total_hits} hits (${ot.high_risk_hits} high-risk)`
      : "No scan";
    setGauge(els.gaugeSpec, g.specificity_score);

    els.valMmej.textContent = `${g.mmej.out_of_frame_probability.toFixed(1)}%`;
    els.gradeMmej.textContent = g.mmej.knockout_grade;
    els.subMmej.textContent = `Bae Score ${g.mmej.microhomology_score.toFixed(0)} · ${g.mmej.total_patterns} MH`;
    setGauge(els.gaugeMmej, g.mmej.out_of_frame_probability);

    // 3. On-Target Biophysical Breakdown
    const od = g.on_target_detail;
    els.onTargetMetaGrid.innerHTML = `
      <span class="meta-chip">GC <b>${od.gc_percent}%</b></span>
      <span class="meta-chip">Tm <b>${od.tm_celsius} °C</b></span>
      <span class="meta-chip">Hairpin <b>${od.hairpin_stem_bp} bp (${od.estimated_hairpin_dg} kcal/mol)</b></span>
      <span class="meta-chip">Scaffold Stem <b>${od.scaffold_stem_bp} bp</b></span>
      <span class="meta-chip">Max Homopolymer <b>${od.homopolymer_run} nt</b></span>
      <span class="meta-chip">Pol III TTTT <b>${od.pol3_terminator ? "YES" : "None"}</b></span>
    `;

    const compNames = I18N[state.lang].compNames;
    els.onTargetComponents.innerHTML = Object.entries(od.components)
      .map(([key, val]) => {
        const pct = Math.min(100, Math.abs(val) * 65);
        const col = val >= 0 ? "#10b981" : "#f43f5e";
        const sign = val > 0 ? "+" : "";
        return `
          <div class="comp-row">
            <div class="comp-top">
              <span>${escapeHtml(compNames[key] || key)}</span>
              <span style="color:${col}">${sign}${val.toFixed(2)}</span>
            </div>
            <div class="comp-track">
              <div class="comp-fill" style="width:${pct}%;background:${col}"></div>
            </div>
          </div>
        `;
      })
      .join("");

    // 4. Interactive Microhomology Deletion Alignment Diagram
    renderMmeJDiagram(g.mmej);

    // 5. Off-Target Mismatch Matrix
    renderOffTargetMatrix(g);
  }

  function renderMmeJDiagram(mmej) {
    els.mmejFsBar.style.width = `${mmej.out_of_frame_probability}%`;
    els.mmejIfBar.style.width = `${mmej.in_frame_probability}%`;
    els.mmejRatioText.textContent = `OOF ${mmej.out_of_frame_probability}% / IF ${mmej.in_frame_probability}%`;

    const patterns = mmej.top_patterns || [];
    if (patterns.length === 0) {
      els.mmejRefTrack.innerHTML = `${escapeHtml(mmej.left_flank)}<span class="mmej-cut-pipe">|</span>${escapeHtml(mmej.right_flank)}`;
      els.mmejPatternsList.innerHTML = `<div class="empty-note">${escapeHtml(t("noMmej"))}</div>`;
      return;
    }

    const activeIdx = Math.min(state.activeMmeJIdx, patterns.length - 1);
    const activePat = patterns[activeIdx];
    const leftLen = mmej.left_flank.length;

    // Highlight BOTH microhomology repeat copies on the reference flanking sequence
    const fullSeq = mmej.left_flank + mmej.right_flank;
    const lStart = activePat.left_mh_start;
    const lEnd = lStart + activePat.mh_length;
    const rStart = activePat.right_mh_start;
    const rEnd = rStart + activePat.mh_length;

    let refHtml = "";
    for (let i = 0; i < fullSeq.length; i++) {
      if (i === leftLen) {
        refHtml += `<span class="mmej-cut-pipe">|</span>`;
      }
      const inLeftMh = i >= lStart && i < lEnd;
      const inRightMh = i >= rStart && i < rEnd;
      if (inLeftMh || inRightMh) {
        refHtml += `<span class="mh-hl">${fullSeq[i]}</span>`;
      } else {
        refHtml += fullSeq[i];
      }
    }
    els.mmejRefTrack.innerHTML = refHtml;

    els.mmejPatternsList.innerHTML = patterns
      .slice(0, 7)
      .map((pat, idx) => {
        const badgeClass = pat.is_frameshift ? "frame-fs" : "frame-if";
        const badgeText = pat.is_frameshift ? t("frameshift") : t("inFrame");
        // Insert visual pipe gap at leftLen to match reference track alignment
        const aLeft = pat.aligned_outcome.slice(0, leftLen);
        const aRight = pat.aligned_outcome.slice(leftLen);
        const colorAligned = (seg) =>
          escapeHtml(seg).replace(/-+/g, (m) => `<span class="mmej-del-dash">${m}</span>`);

        return `
          <div class="mmej-row ${idx === activeIdx ? "active" : ""}" data-mmej-idx="${idx}">
            <div class="mmej-row-head">
              <span>MH <b>${escapeHtml(pat.microhomology)}</b> (${pat.mh_length} bp) · <b>-${pat.deletion_length} bp</b></span>
              <span>
                <span class="frame-badge ${badgeClass}">${escapeHtml(badgeText)}</span>
                <span class="mono-muted" style="margin-left:6px;">Score ${pat.pattern_score.toFixed(0)}</span>
              </span>
            </div>
            <div class="mmej-align-seq">${colorAligned(aLeft)}<span class="mmej-cut-pipe">|</span>${colorAligned(aRight)}</div>
          </div>
        `;
      })
      .join("");

    els.mmejPatternsList.querySelectorAll(".mmej-row").forEach((row) => {
      row.addEventListener("click", () => {
        state.activeMmeJIdx = parseInt(row.getAttribute("data-mmej-idx"), 10);
        renderMmeJDiagram(mmej);
      });
    });
  }

  function renderOffTargetMatrix(g) {
    const ot = g.offtarget_summary;
    if (!ot || !ot.top_hits || ot.top_hits.length === 0) {
      els.offTargetMatrix.innerHTML = `<div class="empty-note">${escapeHtml(t("noOffTargets"))}</div>`;
      return;
    }

    const isCas9 = g.enzyme === "spcas9";
    els.offTargetMatrix.innerHTML = ot.top_hits
      .slice(0, 12)
      .map((hit) => {
        const detailMap = new Map();
        (hit.mismatch_details || []).forEach((d) => detailMap.set(d.position, d));

        const cells = [];
        if (!isCas9) {
          cells.push(
            `<span class="ot-pam ${hit.pam_score < 1 ? "noncanon" : ""}" style="margin-left:0;margin-right:6px;">${escapeHtml(hit.pam)}</span>`
          );
        }
        for (let i = 0; i < hit.alignment.length; i++) {
          const pos = i + 1;
          const ch = hit.alignment[i];
          const inSeed = isCas9 ? pos >= 11 : pos <= 8;
          const det = detailMap.get(pos);
          let cls = "ot-cell";
          let tip = `pos ${pos}: match`;
          if (inSeed) cls += " seed-col";
          if (det) {
            if (det.is_wobble) {
              cls += " mm-wobble";
            } else if (det.is_seed) {
              cls += " mm-seed";
            } else {
              cls += " mm-distal";
            }
            tip = `pos ${pos}: ${det.pair} (${det.is_wobble ? "wobble" : det.is_seed ? "seed" : "distal"}), CFD x${det.cfd_factor}`;
          }
          cells.push(`<span class="${cls}" title="${escapeHtml(tip)}">${escapeHtml(ch)}</span>`);
        }
        if (isCas9) {
          cells.push(
            `<span class="ot-pam ${hit.pam_score < 1 ? "noncanon" : ""}">${escapeHtml(hit.pam)}</span>`
          );
        }

        const riskColor =
          hit.cfd_score >= 0.5
            ? "#f43f5e"
            : hit.cfd_score >= 0.2
            ? "#f59e0b"
            : hit.cfd_score >= 0.05
            ? "#06b6d4"
            : "#94a3b8";

        return `
          <div class="ot-row">
            <div class="ot-row-top">
              <span><b>${escapeHtml(hit.contig)}</b> : ${hit.start + 1}–${hit.end} (${hit.strand})</span>
              <span>
                MM: <b>${hit.mismatches}</b> (Seed ${hit.seed_mismatches}) ·
                CFD: <b style="color:${riskColor}">${hit.cfd_score.toFixed(4)}</b>
              </span>
            </div>
            <div class="ot-tiles">${cells.join("")}</div>
          </div>
        `;
      })
      .join("");
  }

  /* ---------------- Export CSV / FASTA ---------------- */
  function downloadText(filename, content, mime) {
    const blob = new Blob([content], { type: mime });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
  }

  function exportCsv() {
    if (!state.report) return;
    const guides = state.filteredGuides.length ? state.filteredGuides : state.report.guides;
    const header = [
      "rank", "id", "enzyme", "strand", "start_1based", "end_1based", "cut_pos",
      "spacer", "pam", "gc_percent", "on_target_score", "specificity_score",
      "mmej_oof_percent", "composite_score", "warnings",
    ];
    const rows = guides.map((g) =>
      [
        g.rank, g.id, g.enzyme, g.strand, g.start + 1, g.end, g.cut_pos,
        g.spacer, g.pam, g.gc_percent, g.on_target_score, g.specificity_score,
        g.mmej.out_of_frame_probability, g.composite_score, g.warnings.join(";"),
      ].join(",")
    );
    downloadText(
      `crisprguard_${state.report.target_name}_${state.enzyme}.csv`,
      [header.join(","), ...rows].join("\n") + "\n",
      "text/csv;charset=utf-8"
    );
  }

  function exportFasta() {
    if (!state.report) return;
    const guides = state.filteredGuides.length ? state.filteredGuides : state.report.guides;
    const lines = [];
    guides.forEach((g) => {
      lines.push(
        `>${g.id}_${state.report.target_name} enzyme=${g.enzyme} strand=${g.strand} pos=${g.start + 1}-${g.end} cut=${g.cut_pos} pam=${g.pam} composite=${g.composite_score} on_target=${g.on_target_score} specificity=${g.specificity_score} mmej_oof=${g.mmej.out_of_frame_probability}%`
      );
      lines.push(g.spacer);
    });
    downloadText(
      `crisprguard_${state.report.target_name}_${state.enzyme}.fasta`,
      lines.join("\n") + "\n",
      "text/plain;charset=utf-8"
    );
  }

  /* ---------------- Event Binding & Init ---------------- */
  function bindEvents() {
    els.langEnBtn.addEventListener("click", () => {
      state.lang = "en";
      applyLanguage();
    });
    els.langTrBtn.addEventListener("click", () => {
      state.lang = "tr";
      applyLanguage();
    });

    els.enzymeToggle.querySelectorAll(".seg-btn").forEach((btn) => {
      btn.addEventListener("click", () => {
        state.enzyme = btn.getAttribute("data-enzyme");
        els.enzymeToggle.querySelectorAll(".seg-btn").forEach((b) =>
          b.classList.toggle("active", b === btn)
        );
        runDesign();
      });
    });

    els.targetInput.addEventListener("input", () => {
      updateCharCounts();
      if (state.activeSampleId !== "custom") {
        state.activeSampleId = "custom";
        renderPresets();
        updatePresetDescription();
      }
    });
    els.bgInput.addEventListener("input", updateCharCounts);

    els.runBtn.addEventListener("click", runDesign);
    els.exportCsvBtn.addEventListener("click", exportCsv);
    els.exportFastaBtn.addEventListener("click", exportFasta);

    els.zoomSlider.addEventListener("input", () => {
      state.zoom = parseFloat(els.zoomSlider.value);
      els.zoomVal.textContent = `${state.zoom.toFixed(1)}x`;
      drawCanvas();
    });

    els.canvas.addEventListener("mousemove", handleCanvasMove);
    els.canvas.addEventListener("mouseleave", () => {
      state.hoveredGuideId = null;
      els.canvasTooltip.classList.add("hidden");
      drawCanvas();
    });
    els.canvas.addEventListener("click", handleCanvasClick);
    window.addEventListener("resize", () => {
      if (state.report) drawCanvas();
    });

    els.strandFilter.querySelectorAll(".mini-pill").forEach((btn) => {
      btn.addEventListener("click", () => {
        state.strandFilter = btn.getAttribute("data-strand");
        els.strandFilter.querySelectorAll(".mini-pill").forEach((b) =>
          b.classList.toggle("active", b === btn)
        );
        applyFiltersAndRender();
      });
    });

    els.hidePol3Check.addEventListener("change", () => {
      state.hidePol3 = els.hidePol3Check.checked;
      applyFiltersAndRender();
    });

    els.guideSearchInput.addEventListener("input", () => {
      state.searchQuery = els.guideSearchInput.value;
      applyFiltersAndRender();
    });

    els.closeDrawerBtn.addEventListener("click", closeDrawer);
    els.drawerBackdrop.addEventListener("click", closeDrawer);
    window.addEventListener("keydown", (ev) => {
      if (ev.key === "Escape") closeDrawer();
    });

    els.copySpacerBtn.addEventListener("click", async () => {
      if (!state.report || !state.selectedGuideId) return;
      const g = state.report.guides.find((x) => x.id === state.selectedGuideId);
      if (!g) return;
      try {
        await navigator.clipboard.writeText(g.spacer);
        els.copySpacerBtn.textContent = t("copied");
        setTimeout(() => {
          els.copySpacerBtn.textContent = t("copySpacer");
        }, 1200);
      } catch {
        // Clipboard unavailable in insecure context
      }
    });
  }

  async function init() {
    bindEvents();
    try {
      const resp = await fetch("/api/samples");
      const data = await resp.json();
      state.samples = data.samples || [];
      applyLanguage();
      if (state.samples.length > 0) {
        selectSample(state.samples[0].id, true);
      }
    } catch (err) {
      els.statusMessage.textContent = `Failed to load samples: ${err.message}`;
      els.statusMessage.classList.add("error");
    }
  }

  init();
})();
