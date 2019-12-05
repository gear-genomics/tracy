module.exports = {
  title: "Tracy documentation",
  description:
    "Documentation of Tracy, an app for basecalling, alignment, assembly and deconvolution of Sanger Chromatogram trace files",
  base: "/docs/tracy/",
  themeConfig: {
    repo: "gear-genomics/tracy",
    nav: [
      { text: "Home", link: "/" },
      { text: "Installation", link: "/installation/" },
      { text: "Usage", link: "/cli/" },
      { text: "Web Apps", link: "/webapps/" },
      { text: "FAQ", link: "/faq/" }
    ],
    sidebar: ["/installation/", "/cli/", "/webapps/", "/faq/"]
  },
  plugins: {
    "@vuepress/back-to-top": true
  }
};
