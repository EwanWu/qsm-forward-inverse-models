'use strict';
const fs = require('node:fs');
const path = require('node:path');
const Q = require('../src/qsm-core.js');

const result = Q.run('forward');
const source = result.maps.find(x => x.key === 'susceptibility').data;
const field = result.maps.find(x => x.key === 'localField').data;
function panel(values, x0, label) {
  const max = Math.max(...values.map(Math.abs)) || 1;
  const cells = Array.from(values, (v, i) => {
    const x = i % 32, y = Math.floor(i / 32), a = Math.abs(v) / max;
    const color = v >= 0 ? `rgb(${Math.round(255 - 105 * a)},${Math.round(255 - 175 * a)},${Math.round(255 - 195 * a)})`
                         : `rgb(${Math.round(255 - 195 * a)},${Math.round(255 - 170 * a)},${Math.round(255 - 90 * a)})`;
    return `<rect x="${x0 + x * 8}" y="${45 + y * 8}" width="8" height="8" fill="${color}"/>`;
  }).join('');
  return `<text x="${x0}" y="30" font-family="Arial" font-size="16">${label}</text>${cells}`;
}
const svg = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 560 332" role="img" aria-label="Synthetic susceptibility source and dipole field"><rect width="560" height="332" fill="white"/>${panel(source, 16, 'Synthetic susceptibility')}${panel(field, 288, 'Dipole field')}</svg>`;
const output = path.join(__dirname, 'qsm_forward_example.svg');
fs.writeFileSync(output, svg);
console.log(JSON.stringify({grid: '32 x 32 x 32', centerField: result.stats.centerField, figure: output}, null, 2));
