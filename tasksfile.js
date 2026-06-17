const { sh, cli } = require("tasksfile");

function cleanAll() {
    sh("rm -rf build");
}

function downloadGoogleTest() {
    sh("mkdir -p build");
    sh("wget https://github.com/google/googletest/archive/release-1.10.0.tar.gz", { cwd: "build" });
    sh("tar xzf release-1.10.0.tar.gz", { cwd: "build" });
    sh("rm release-1.10.0.tar.gz", { cwd: "build" });
}

function compileGoogleTest() {
    sh("g++ -Igoogletest -Igoogletest/include -c googletest/src/gtest-all.cc", {
        cwd: "build/googletest-release-1.10.0",
        nopipe: true
    });
    sh("ar -rv libgtest.a gtest-all.o", {
        cwd: "build/googletest-release-1.10.0",
        nopipe: true
    });
}

function createMpSources() {
    sh("mkdir -p build");

    const renderMpHpp =
        String.raw`node -e "const fs=require('fs'); const ejs=require('ejs'); const path=require('path'); const data={ n64: 4 }; const tpl=fs.readFileSync(path.join('..','src','mp.hpp.ejs'),'utf8'); const out=ejs.render(tpl, data); fs.writeFileSync('mp.hpp', out);"`;

    const renderMpCpp =
        String.raw`node -e "const fs=require('fs'); const ejs=require('ejs'); const path=require('path'); const data={ n64: 4 }; const tpl=fs.readFileSync(path.join('..','src','mp.cpp.ejs'),'utf8'); const out=ejs.render(tpl, data); fs.writeFileSync('mp.cpp', out);"`;

    sh(renderMpHpp, { cwd: "build", nopipe: true });
    sh(renderMpCpp, { cwd: "build", nopipe: true });
}

function createFieldSources() {
    sh("mkdir -p build");

    sh(
        "node ../src/buildzqfield.js -q 21888242871839275222246405745257275088696311157297823662689037894645226208583 -n Fq",
        { cwd: "build", nopipe: true }
    );
    sh(
        "node ../src/buildzqfield.js -q 21888242871839275222246405745257275088548364400416034343698204186575808495617 -n Fr",
        { cwd: "build", nopipe: true }
    );

    if (process.platform === "linux" && process.arch === "x64") {
        sh("nasm -felf64 fq.asm -o fq_asm.o", { cwd: "build", nopipe: true });
        sh("nasm -felf64 fr.asm -o fr_asm.o", { cwd: "build", nopipe: true });
    } else if (process.platform === "darwin" && process.arch === "x64") {
        sh("nasm -fmacho64 --prefix _ fq.asm -o fq_asm.o", { cwd: "build", nopipe: true });
        sh("nasm -fmacho64 --prefix _ fr.asm -o fr_asm.o", { cwd: "build", nopipe: true });
    }
}

function getIncludeFlags(extra = []) {
    const flags = [
        "-Igoogletest-release-1.10.0/googletest/include",
        "-I.",
        ...extra
    ];
    return flags.join(" ");
}

function getFieldSourcesForTest() {
    if (process.platform === "darwin" && process.arch === "arm64") {
        return [
            "fq.cpp",
            "fq_generic.cpp",
            "fq_raw_arm64.s",
            //"fq_raw_generic.cpp",
            "fr.cpp",
            "fr_generic.cpp",
            "fr_raw_arm64.s"
            //"fr_raw_generic.cpp"
        ];
    }

    if (process.arch === "x64") {
        const asmExists = process.platform === "linux" || process.platform === "darwin";
        if (asmExists) {
            return [
                "fq.cpp",
                "fq_asm.o",
                "fr.cpp",
                "fr_asm.o"
            ];
        }
    }

    return [
        "fq.cpp",
        "fq_generic.cpp",
        "fq_raw_generic.cpp",
        "fr.cpp",
        "fr_generic.cpp",
        "fr_raw_generic.cpp"
    ];
}

function getFieldSourcesForBench() {
    if (process.platform === "darwin" && process.arch === "arm64") {
        return [
            "fq.cpp",
            "fq_generic.cpp",
            "fq_raw_arm64.s",
            "fr.cpp",
            "fr_generic.cpp",
            "fr_raw_arm64.s"
        ];
    }

    if (process.arch === "x64") {
        const asmExists = process.platform === "linux" || process.platform === "darwin";
        if (asmExists) {
            return [
                "fq.cpp",
                "fq_asm.o",
                "fr.cpp",
                "fr_asm.o"
            ];
        }
    }

    return [
        "fq.cpp",
        "fq_generic.cpp",
        "fq_raw_generic.cpp",
        "fr.cpp",
        "fr_generic.cpp",
        "fr_raw_generic.cpp"
    ];
}

function getOpenMpFlags() {
    if (process.platform === "linux") {
        return "-fopenmp";
    }
    return "";
}

function getAsmDefines() {
    if (process.platform === "darwin" && process.arch === "arm64") {
        return "-DUSE_ASM -DARCH_ARM64";
    }
    if (process.arch === "x64" && (process.platform === "darwin" || process.platform === "linux")) {
        return "-DUSE_ASM -DARCH_X86_64";
    }
    return "";
}

function testSplitParStr() {
    const asmDefines = getAsmDefines();

    const cmd =
        "g++ " +
        (asmDefines ? asmDefines + " " : "") +
        getIncludeFlags(["-I../src"]) +
        " ../c/splitparstr.cpp" +
        " ../c/splitparstr_test.cpp" +
        " googletest-release-1.10.0/libgtest.a" +
        " -pthread -std=c++11 -o splitparsestr_test";

    sh(cmd, { cwd: "build", nopipe: true });
    sh("./splitparsestr_test", { cwd: "build", nopipe: true });
}

function testAltBn128() {
    const fieldSources = getFieldSourcesForTest().join(" ");
    const asmDefines = getAsmDefines();

    const cmd =
        "g++ -O3 " +
        (asmDefines ? asmDefines + " " : "") +
        getIncludeFlags(["-I../c"]) +
        " ../c/naf.cpp" +
        " ../c/splitparstr.cpp" +
        " ../c/alt_bn128.cpp" +
        " ../c/alt_bn128_test.cpp" +
        " ../c/misc.cpp" +
        " mp.cpp" +
        " " + fieldSources +
        " googletest-release-1.10.0/libgtest.a" +
        " -o altbn128_test" +
        " -pthread -std=c++17 -g";

    sh(cmd, { cwd: "build", nopipe: true });
    sh("./altbn128_test", { cwd: "build", nopipe: true });
}

function benchMultiExpG1() {
    const fieldSources = getFieldSourcesForBench().join(" ");
    const omp = getOpenMpFlags();
    const asmDefines = getAsmDefines();

    const cmd =
        "g++ -O3 -g " +
        (asmDefines ? asmDefines + " " : "") +
        getIncludeFlags(["-I../c"]) +
        " ../c/naf.cpp" +
        " ../c/splitparstr.cpp" +
        " ../c/alt_bn128.cpp" +
        " ../c/misc.cpp" +
        " ../benchmark/multiexp_g1.cpp" +
        " mp.cpp" +
        " " + fieldSources +
        " -o multiexp_g1_benchmark " +
        " -pthread -std=c++11" +
        (omp ? " " + omp : "");

    sh(cmd, { cwd: "build", nopipe: true });
    sh("./multiexp_g1_benchmark 16777216", { cwd: "build", nopipe: true });
}

function benchMultiExpG2() {
    const fieldSources = getFieldSourcesForBench().join(" ");
    const omp = getOpenMpFlags();
    const asmDefines = getAsmDefines();

    const cmd =
        "g++ -O3 -g " +
        (asmDefines ? asmDefines + " " : "") +
        getIncludeFlags(["-I../c"]) +
        " ../c/naf.cpp" +
        " ../c/splitparstr.cpp" +
        " ../c/alt_bn128.cpp" +
        " ../c/misc.cpp" +
        " ../benchmark/multiexp_g2.cpp" +
        " mp.cpp" +
        " " + fieldSources +
        " -o multiexp_g2_benchmark " +
        " -pthread -std=c++11" +
        (omp ? " " + omp : "");

    sh(cmd, { cwd: "build", nopipe: true });
    sh("./multiexp_g2_benchmark 16777216", { cwd: "build", nopipe: true });
}

cli({
    cleanAll,
    downloadGoogleTest,
    compileGoogleTest,
    createMpSources,
    createFieldSources,
    testSplitParStr,
    testAltBn128,
    benchMultiExpG1,
    benchMultiExpG2,
});
