/* Data entry point for the shared proton tag-and-probe processor. */
String[] forwarded = args
GroovyClassLoader loader = new GroovyClassLoader(this.class.classLoader)
loader.parseClass(new File('processing_scripts/proton_efficiency_common.groovy'))
loader.loadClass('ProtonEfficiencyProcessor').run(forwarded, false)

