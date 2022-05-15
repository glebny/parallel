import qbs

CppApplication {
    consoleApplication: true
    install: true
    files: [
        "blood.c",
        "blood.h",
        "brusselator.c",
        "brusselator.h",
        "common.c",
        "common.h",
        "cursedpart.c",
        "cursedpart.h",
        "main.c",
        "reduction.c",
        "reduction.h",
        "vtk.c",
        "vtk.h",
    ]
}
