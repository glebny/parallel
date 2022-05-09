import qbs

CppApplication {
    consoleApplication: true
    install: true
    files: [
        "brusselator.h",
        "common.h",
        "main.c",
        "reduction.h",
        "vtk.h",
    ]
}
