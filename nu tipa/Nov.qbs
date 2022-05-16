import qbs

CppApplication {
    consoleApplication: true
    files: [
        "../parallel-main/NOv/README.txt",
        "../parallel-main/NOv/SRC.qbs",
        "../parallel-main/NOv/SRC.qbs.user",
        "../parallel-main/NOv/blood.c",
        "../parallel-main/NOv/blood.h",
        "../parallel-main/NOv/brusselator.c",
        "../parallel-main/NOv/brusselator.h",
        "../parallel-main/NOv/common.c",
        "../parallel-main/NOv/common.h",
        "../parallel-main/NOv/cursedpart.c",
        "../parallel-main/NOv/cursedpart.h",
        "../parallel-main/NOv/curses.h",
        "../parallel-main/NOv/main.c",
        "../parallel-main/NOv/pdcurses.a",
        "../parallel-main/NOv/reduction.c",
        "../parallel-main/NOv/reduction.h",
        "../parallel-main/NOv/vtk.c",
        "../parallel-main/NOv/vtk.h",
        "blood.c",
        "blood.h",
        "brusselator.c",
        "brusselator.h",
        "common.c",
        "common.h",
        "cursedpart.c",
        "cursedpart.h",
        "curses.h",
        "main.c",
        "pdcurses.a",
        "reduction.c",
        "reduction.h",
        "vtk.c",
        "vtk.h",
    ]

    Group {     // Properties for the produced executable
        fileTagsFilter: "application"
        qbs.install: true
        qbs.installDir: "bin"
    }
}
