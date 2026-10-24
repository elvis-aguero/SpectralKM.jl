#!/usr/bin/env python3
"""
Audits assets/schematic.svg for overlaps that should not exist, by measuring them in a real browser layout.

Every <text> label is tested against (i) every other label (bounding boxes) and (ii) every stroked element
(lines, paths, circles) by sampling a grid of points inside the label's bounding box with isPointInStroke.
A label sitting on a fill (inside the drop, inside the bath) is fine and is not flagged. Exit status 1 on any
collision.

    python3 derivations/audit_schematic.py [path/to/schematic.svg]

Leader lines are allowed to END at a label if they are listed in ALLOW (substring of the label, substring of the
element's `data-note`), but the generator does not use that escape hatch: it leaves a gap instead.
"""
import json, pathlib, re, subprocess, sys, tempfile

CHROME = "/opt/pw-browsers/chromium-1243/chrome-linux64/chrome"
svg_path = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else pathlib.Path(__file__).resolve().parent.parent / "assets" / "schematic.svg"
svg = svg_path.read_text()

js = r"""
const svg = document.querySelector('svg');
const texts = [...svg.querySelectorAll('text')];
const shapes = [...svg.querySelectorAll('line, path, circle')].filter(e => {
  const s = e.getAttribute('stroke'); return s && s !== 'none' && e.getAttribute('stroke') !== '#fff';
});
const out = [];
const boxes = texts.map(t => { const b = t.getBBox(); return {t, x: b.x, y: b.y + 0.18*b.height, w: b.width, h: 0.66*b.height, s: t.textContent}; });
for (let i = 0; i < boxes.length; i++) {
  const a = boxes[i];
  for (let j = i + 1; j < boxes.length; j++) {
    const b = boxes[j];
    if (a.x < b.x + b.w && b.x < a.x + a.w && a.y < b.y + b.h && b.y < a.y + a.h)
      out.push('TEXT/TEXT   "' + a.s + '"  x  "' + b.s + '"');
  }
  const pts = [];
  for (let u = 0; u <= 16; u++) for (let v = 0; v <= 4; v++) pts.push([a.x + a.w * u / 16, a.y + a.h * v / 4]);
  for (const e of shapes) {
    const sw = parseFloat(e.getAttribute('stroke-width') || '1');
    e.style.strokeWidth = sw;  // isPointInStroke uses the element's stroke-width
    let hit = 0;
    for (const [x, y] of pts) { const p = new DOMPoint(x, y); if (e.isPointInStroke(p)) hit++; }
    if (hit > 0) {
      const d = e.tagName === 'line' ? `line (${e.getAttribute('x1')},${e.getAttribute('y1')})->(${e.getAttribute('x2')},${e.getAttribute('y2')})`
              : e.tagName === 'circle' ? `circle r=${e.getAttribute('r')}` : 'path ' + (e.getAttribute('d') || '').slice(0, 40);
      out.push('TEXT/STROKE "' + a.s + '"  x  ' + d + '  (' + hit + ' of ' + pts.length + ' sample points)');
    }
  }
}
document.getElementById('audit_out').textContent = JSON.stringify(out);
"""
html = f"<!doctype html><html><body style='margin:0'>{svg}<pre id='audit_out'></pre><script>window.onload=function(){{{js}}}</script></body></html>"
with tempfile.TemporaryDirectory() as d:
    f = pathlib.Path(d) / "a.html"; f.write_text(html)
    res = subprocess.run([CHROME, "--headless", "--no-sandbox", "--disable-gpu", "--virtual-time-budget=4000", "--dump-dom", f"file://{f}"],
                         capture_output=True, text=True).stdout
m = re.search(r'<pre id="audit_out">(.*?)</pre>', res, re.S)
if not m:
    print("could not read the audit result from the browser"); sys.exit(2)
issues = json.loads(m.group(1).replace("&quot;", '"').replace("&amp;", "&").replace("&gt;", ">").replace("&lt;", "<"))
for i in issues:
    print(i)
print(f"{len(issues)} collision(s) in {svg_path.name}")
sys.exit(1 if issues else 0)
