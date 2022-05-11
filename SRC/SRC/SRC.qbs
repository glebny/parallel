import qbs

CppApplication {
    consoleApplication: true
    install: true
    files: [
        "brusselator.h",
        "common.h",
        "cursedpart.c",
        "cursedpart.h",
        "main.c",
        "reduction.c",
        "reduction.h",
        "vtk.h",
    ]
}
