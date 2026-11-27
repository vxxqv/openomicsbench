# Website views

GitHub repository Traffic measures visits to the GitHub repository. It does not count visits to the [OpenOmicsBench website](https://vxxqv.github.io/openomicsbench/).

The GitHub Pages deployment can optionally add [Cloudflare Web Analytics](https://developers.cloudflare.com/web-analytics/get-started/) to its HTML pages. This is a free, privacy-focused page-view service that works on GitHub Pages without moving the site's DNS to Cloudflare. The tracking code is added only when the GitHub Actions repository variable `CLOUDFLARE_WEB_ANALYTICS_TOKEN` is set; an empty variable leaves the website unchanged.

## One-time setup

1. Sign in to a Cloudflare account and open **Web Analytics**.
2. Select **Add a site** and enter `vxxqv.github.io` as the hostname. The site path is `/openomicsbench/`; the manual beacon is installed only on this project's pages.
3. Copy the JavaScript snippet from **Manage site**. Its `data-cf-beacon` attribute contains a 32-character site token.
4. In the GitHub repository, open **Settings → Secrets and variables → Actions → Variables → New repository variable**. Name it `CLOUDFLARE_WEB_ANALYTICS_TOKEN` and paste only the token as the value.
5. After the analytics change reaches `main`, the **Publish project site** workflow deploys automatically. Confirm the deployed homepage contains the Cloudflare beacon and that a real visit appears in the Cloudflare dashboard after a few minutes. Any additional HTML pages assembled by this workflow also receive the beacon.

The token is a public site identifier embedded in page HTML, not an API password. The build script checks its format and refuses a malformed value. Analytics starts collecting after the beacon is deployed; it cannot reconstruct earlier website visits. Cloudflare's dashboard reports page views and visitors; neither figure is a count of GitHub clones, package downloads or distinct users of the software.
