import { defineUserConfig } from '@vuepress/cli'
import { viteBundler } from '@vuepress/bundler-vite'
import { defaultTheme } from '@vuepress/theme-default'

export default defineUserConfig({
  title: "Tracy documentation",
  description: "Documentation of Tracy, an app for basecalling, alignment, assembly and deconvolution of Sanger Chromatogram trace files",
  base: "/docs/tracy/",
  theme: defaultTheme({
    repo: "gear-genomics/tracy",
    navbar: [
      { text: "Home", link: "/" },
      { text: "Installation", link: "/installation/" },
      { text: "Usage", link: "/cli/" },
      { text: "Web Apps", link: "/webapp/" },
      { text: "FAQ", link: "/faq/" }
    ],
    sidebar: false
  }),
  bundler: viteBundler({})
})
