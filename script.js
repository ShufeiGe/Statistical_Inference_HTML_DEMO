function initDemoNavigation() {
  const links = Array.from(document.querySelectorAll("a[href]"));

  // 没有导航链接时，不添加搜索框
  if (links.length === 0) return;

  const searchInput = document.createElement("input");
  searchInput.id = "demo-search";
  searchInput.type = "search";
  searchInput.placeholder = "搜索 demo，例如：置信区间";
  searchInput.setAttribute("aria-label", "搜索演示");

  const status = document.createElement("p");
  status.id = "search-status";
  status.setAttribute("role", "status");

  // 在第一个链接前插入搜索框
  links[0].before(searchInput);

  // 在最后一个链接后插入筛选结果提示
  links[links.length - 1].after(status);

  function filterDemos() {
    const keyword = searchInput.value.trim().toLowerCase();
    let visibleCount = 0;

    links.forEach(function (link) {
      const title = link.textContent.toLowerCase();
      const matches = title.includes(keyword);

      link.hidden = !matches;

      if (matches) {
        visibleCount++;
      }
    });

    if (visibleCount === 0) {
      status.textContent = "没有找到匹配的 demo，请尝试其他关键词。";
    } else {
      status.textContent = `显示 ${visibleCount} / ${links.length} 个 demo`;
    }
  }

  searchInput.addEventListener("input", filterDemos);
  filterDemos();
}

if (document.readyState === "loading") {
  document.addEventListener("DOMContentLoaded", initDemoNavigation);
} else {
  initDemoNavigation();
}
