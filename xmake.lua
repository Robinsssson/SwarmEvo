add_rules("mode.debug", "mode.release")
add_rules("plugin.compile_commands.autoupdate", {outputdir = "."})

includes("algmath")

set_toolchains("gcc")

target("alg-project")
    set_kind("shared")
    add_files("src/**.c")
    add_defines("ALG_EXPORT")
    add_deps("algmath")
    add_includedirs("algmath")
    add_cflags(
        "-Wall", "-Wextra", "-Werror", "-pedantic", "-std=c11",
        "-Wshadow", "-Wconversion", "-Wfloat-equal", "-Wundef", 
        "-Wstrict-prototypes", "-Wmissing-prototypes", 
        "-Wredundant-decls"
    )

target("test-alg-project")
    set_kind("binary")
    add_files("test/**.c")
    add_deps("alg-project")
    add_cflags(
        "-Wall", "-Wextra", "-Werror", "-pedantic", "-std=c11",
        "-Wshadow", "-Wconversion", "-Wfloat-equal", "-Wundef", 
        "-Wstrict-prototypes", "-Wmissing-prototypes", 
        "-Wredundant-decls"
    )
