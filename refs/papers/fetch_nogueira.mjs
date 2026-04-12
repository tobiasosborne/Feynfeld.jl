#!/usr/bin/env node
// Fetch Nogueira 1993 (J. Comput. Phys. 105, 279) via Playwright + TIB VPN.
// Reuses the qvls-sturm playwright install and Sturm.jl browser profile.
//
// IMPORTANT: this script waits up to 5 minutes for you to:
//   - solve any CAPTCHA
//   - accept cookie/consent dialogs
//   - wait for TIB institutional access to register
// It polls for citation_pdf_url meta tag or a Download PDF button
// as the signal that institutional access is live.

import { chromium } from '/home/tobiasosborne/Projects/qvls-sturm/viz/node_modules/playwright/index.mjs';
import { writeFileSync, existsSync, mkdirSync } from 'fs';
import { resolve } from 'path';

const BASE = resolve(import.meta.dirname);
const OUT = resolve(BASE, 'Nogueira1993_JCompPhys105_279.pdf');
const TRIGGER_URL = 'https://www.sciencedirect.com/science/article/pii/S0021999183710740';
const PDF_URL = 'https://www.sciencedirect.com/science/article/pii/S0021999183710740/pdfft';
const WAIT_MS = 300000; // 5 min for you to solve captcha/login

async function waitForRealAccess(page, timeoutMs) {
  console.log(`\n>>> YOU HAVE ${timeoutMs / 1000}s to:`);
  console.log('    - solve any CAPTCHA in the browser window');
  console.log('    - accept cookie / consent dialogs');
  console.log('    - confirm you see the ScienceDirect article page with a Download PDF button');
  console.log('    The script will auto-detect and proceed.\n');

  try {
    await page.waitForFunction(
      () => {
        const meta = document.querySelector('meta[name="citation_pdf_url"]');
        if (meta && meta.content) return true;
        const btn = Array.from(document.querySelectorAll('a, button')).find((el) =>
          /download\s*pdf|view\s*pdf/i.test(el.textContent || ''),
        );
        if (btn) return true;
        return false;
      },
      { timeout: timeoutMs },
    );
    return true;
  } catch (_) {
    return false;
  }
}

async function main() {
  if (existsSync(OUT)) {
    console.log(`SKIP: already have ${OUT}`);
    return;
  }
  const userDataDir = resolve('/home/tobiasosborne/Projects/Sturm.jl/docs/literature/quantum_simulation', '.browser-profile');
  mkdirSync(userDataDir, { recursive: true });

  console.log('Launching HEADED Chromium...');
  const context = await chromium.launchPersistentContext(userDataDir, {
    headless: false,
    args: ['--disable-blink-features=AutomationControlled'],
    viewport: { width: 1280, height: 900 },
  });
  const page = context.pages()[0] || (await context.newPage());

  console.log(`Trigger: ${TRIGGER_URL}`);
  await page.goto(TRIGGER_URL, { waitUntil: 'domcontentloaded', timeout: 60000 });

  const ok = await waitForRealAccess(page, WAIT_MS);
  if (!ok) {
    console.log('TIMEOUT: institutional access not confirmed within 5 minutes.');
    console.log('Leaving browser open. You can still navigate manually.');
    return;
  }
  console.log('Access confirmed.');
  await new Promise((r) => setTimeout(r, 2000));

  const info = await page.evaluate(() => ({
    title: document.title,
    url: location.href,
    citationPdf: document.querySelector('meta[name="citation_pdf_url"]')?.content,
    pdfLinks: Array.from(document.querySelectorAll('a[href]'))
      .map((a) => a.href)
      .filter((h) => h.includes('pdfft') || h.includes('/pdf/'))
      .slice(0, 10),
  }));
  console.log('Info:', JSON.stringify(info, null, 2));

  let body = null;

  // Strategy 1: citation_pdf_url
  if (info.citationPdf) {
    console.log(`Try citation_pdf_url: ${info.citationPdf}`);
    const r = await page.request.get(info.citationPdf, {
      timeout: 60000,
      headers: { Referer: TRIGGER_URL },
    });
    console.log(`  status=${r.status()}`);
    if (r.status() === 200) {
      const b = await r.body();
      if (b.slice(0, 5).toString() === '%PDF-') body = b;
    }
  }

  // Strategy 2: discovered PDF links
  if (!body) {
    for (const href of info.pdfLinks) {
      console.log(`Try link: ${href}`);
      const r = await page.request.get(href, {
        timeout: 60000,
        headers: { Referer: TRIGGER_URL },
      });
      console.log(`  status=${r.status()}`);
      if (r.status() === 200) {
        const b = await r.body();
        if (b.slice(0, 5).toString() === '%PDF-') {
          body = b;
          break;
        }
      }
    }
  }

  // Strategy 3: direct pdfft endpoint
  if (!body) {
    console.log(`Try direct: ${PDF_URL}`);
    const r = await page.request.get(PDF_URL, {
      timeout: 60000,
      headers: { Referer: TRIGGER_URL },
    });
    console.log(`  status=${r.status()}`);
    if (r.status() === 200) {
      const b = await r.body();
      if (b.slice(0, 5).toString() === '%PDF-') body = b;
    }
  }

  // Strategy 4: click Download PDF button and capture download
  if (!body) {
    console.log('Trying click Download PDF button');
    try {
      const [download] = await Promise.all([
        page.waitForEvent('download', { timeout: 60000 }).catch(() => null),
        page
          .click('a:has-text("Download PDF"), button:has-text("Download PDF"), a:has-text("View PDF")')
          .catch(() => null),
      ]);
      if (download) {
        const p = await download.path();
        if (p) {
          const fs = await import('fs');
          const b = fs.readFileSync(p);
          if (b.slice(0, 5).toString() === '%PDF-') body = b;
        }
      }
    } catch (e) {
      console.log(`  click error: ${e.message}`);
    }
  }

  if (body) {
    writeFileSync(OUT, body);
    console.log(`OK: ${body.length} bytes → ${OUT}`);
  } else {
    console.log('FAIL: could not obtain PDF. Browser left open for manual rescue.');
    // do NOT close context so the user can grab the PDF manually
    return;
  }

  await context.close();
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
