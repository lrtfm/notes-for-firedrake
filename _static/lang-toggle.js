/* Language toggle for bilingual content.
 *
 * Convention:
 *   - blocks with class `lang-zh` are shown only in Chinese mode;
 *   - blocks with class `lang-en` are shown only in English mode;
 *   - unmarked content (all code cells, untranslated text) is always shown.
 *
 * In MyST markdown, mark a block with a colon fence:
 *   :::lang-en
 *   English text ...
 *   :::
 * or add the class to an admonition: `:class: tip lang-en`.
 *
 * The preference is stored in localStorage and applied before first paint
 * to avoid flashing.
 */
(function () {
  var KEY = "firedrake-notes-lang";
  var lang = "zh";
  try {
    lang = localStorage.getItem(KEY) || "zh";
  } catch (e) {
    /* localStorage unavailable (e.g. file://) - keep default */
  }
  if (lang !== "zh" && lang !== "en") {
    lang = "zh";
  }
  var root = document.documentElement;
  root.setAttribute("data-content-lang", lang);

  /* Sidebar part captions come from _toc.yml as plain text, so they are
   * translated here instead of with CSS classes. */
  var CAPTIONS = {
    "Firedrake 入门 (Tutorials)": { zh: "Firedrake 入门", en: "Tutorials" },
    "进阶专题 (Advanced Topics)": { zh: "进阶专题", en: "Advanced Topics" },
    "内部机制 (Internals)": { zh: "内部机制", en: "Internals" },
    "Installation": { zh: "安装", en: "Installation" },
    "Additional Information": { zh: "附加信息", en: "Additional Information" }
  };

  function updateCaptions() {
    var nodes = document.querySelectorAll(
      ".bd-sidebar-primary p.caption > .caption-text, " +
        ".bd-sidebar-primary p.caption"
    );
    nodes.forEach(function (el) {
      if (el.children.length > 0) return; /* handle only text-bearing node */
      var orig = el.getAttribute("data-caption-orig") || el.textContent.trim();
      el.setAttribute("data-caption-orig", orig);
      var tr = CAPTIONS[orig];
      if (tr) {
        el.textContent = tr[lang] || orig;
      }
    });
  }

  function updateDocTitle() {
    /* The <title> tag concatenates both languages at build time; rebuild it
     * from the visible heading text. innerText ignores display:none spans,
     * but the print-only copy of the title is plain text, so skip it. */
    var h1 = Array.prototype.find.call(
      document.querySelectorAll("h1"),
      function (el) {
        return !el.closest("#jb-print-docs-body") && !el.closest(".onlyprint");
      }
    );
    if (!h1) return;
    var text = h1.innerText
      .replace(/#\s*$/, "") /* headerlink anchor */
      .replace(/^\d+(\.\d+)*\.\s*/, "") /* section number */
      .trim();
    if (!text) return;
    document.title =
      text === "Firedrake Notes" ? text : text + " — Firedrake Notes";
  }

  function apply(next) {
    lang = next;
    root.setAttribute("data-content-lang", lang);
    try {
      localStorage.setItem(KEY, lang);
    } catch (e) {}
    var btn = document.getElementById("lang-toggle");
    if (btn) {
      /* The label shows the language you will switch TO. */
      btn.textContent = lang === "zh" ? "EN" : "中文";
      btn.setAttribute(
        "title",
        lang === "zh" ? "Switch to English" : "切换到中文"
      );
      btn.setAttribute(
        "aria-label",
        lang === "zh" ? "Switch to English" : "切换到中文"
      );
    }
    updateCaptions();
    updateDocTitle();
  }

  document.addEventListener("DOMContentLoaded", function () {
    var btn = document.createElement("button");
    btn.id = "lang-toggle";
    btn.type = "button";
    btn.addEventListener("click", function () {
      apply(lang === "zh" ? "en" : "zh");
    });
    var host = document.querySelector(".article-header-buttons");
    if (host) {
      btn.className = "btn btn-sm lang-toggle-btn";
      host.insertBefore(btn, host.firstChild);
    } else {
      btn.className = "lang-toggle-btn lang-toggle-floating";
      document.body.appendChild(btn);
    }
    apply(lang);
  });
})();
