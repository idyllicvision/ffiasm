const chai = require("chai");
const assert = chai.assert;

const fs = require("fs");
const tmp = require("tmp-promise");
const path = require("path");
const util = require("util");
const ejs = require("ejs");
const exec = util.promisify(require("child_process").exec);

const buildZqField = require("../../index.js").buildZqField;

module.exports = testField;

async function renderTemplate(srcPath, data) {
    const tpl = await fs.promises.readFile(srcPath, "utf8");
    return ejs.render(tpl, data);
}

async function testField(prime, test, options = {}) {
    tmp.setGracefulCleanup();

    const dir = await tmp.dir({ prefix: "ffiasm_", unsafeCleanup: true });
    const useAsm = options.useAsm !== false;

    const source = await buildZqField(prime, "Fr");

    await fs.promises.writeFile(path.join(dir.path, "fr.asm"), source.asm, "utf8");
    await fs.promises.writeFile(path.join(dir.path, "fr.hpp"), source.hpp, "utf8");
    await fs.promises.writeFile(path.join(dir.path, "fr.cpp"), source.cpp, "utf8");
    await fs.promises.writeFile(path.join(dir.path, "fr_element.hpp"), source.element_hpp, "utf8");
    await fs.promises.writeFile(path.join(dir.path, "fr_generic.cpp"), source.generic_cpp, "utf8");
    await fs.promises.writeFile(path.join(dir.path, "fr_raw_generic.cpp"), source.raw_generic_cpp, "utf8");
    await fs.promises.writeFile(path.join(dir.path, "fr_raw_arm64.s"), source.raw_arm64_s, "utf8");

    const mpData = { n64: source.n64 };
    //console.log("prime bits =", prime.bitLength().toString());
    //console.log("source.n64 =", source.n64);
    //console.log("mpData =", mpData);
    const mpHpp = await renderTemplate(path.join(__dirname, "..", "..", "src", "mp.hpp.ejs"), mpData);
    //console.log(mpHpp.match(/#define MP_N64 .*/)[0]);
    const mpCpp = await renderTemplate(path.join(__dirname, "..", "..", "src", "mp.cpp.ejs"), mpData);

    await fs.promises.writeFile(path.join(dir.path, "mp.hpp"), mpHpp, "utf8");
    await fs.promises.writeFile(path.join(dir.path, "mp.cpp"), mpCpp, "utf8");

    await exec(`cp ${path.join(__dirname, "tester.cpp")} ${dir.path}`);

    let defines = [];
    let sources = [
        path.join(dir.path, "tester.cpp"),
        path.join(dir.path, "mp.cpp"),
        path.join(dir.path, "fr.cpp")
    ];

    /*if (useAsm && process.platform === "darwin" && process.arch === "arm64") {
        defines.push("-DUSE_ASM", "-DARCH_ARM64");
        sources.push(
            path.join(dir.path, "fr_generic.cpp"),
            path.join(dir.path, "fr_raw_arm64.s")
        );
    } else */{
        sources.push(
            path.join(dir.path, "fr_generic.cpp"),
            path.join(dir.path, "fr_raw_generic.cpp")
        );
    }

    const compileCmd = [
        "g++",
        ...defines,
        `-I${dir.path}`,
        ...sources,
        "-o", path.join(dir.path, "tester"),
        "-std=c++17",
        "-g"
    ].join(" ");

    await exec(compileCmd);

    const inLines = [];
    for (let i = 0; i < test.length; i++) {
        for (let j = 0; j < test[i][0].length; j++) {
            inLines.push(test[i][0][j]);
        }
    }
    inLines.push("");

    await fs.promises.writeFile(path.join(dir.path, "in.tst"), inLines.join("\n"), "utf8");

    await exec(
        `${path.join(dir.path, "tester")} <${path.join(dir.path, "in.tst")} >${path.join(dir.path, "out.tst")}`
    );

    const res = await fs.promises.readFile(path.join(dir.path, "out.tst"), "utf8");
    const resLines = res.split("\n");

    for (let i = 0; i < test.length; i++) {
        const expected = test[i][1].toString();
        const calculated = resLines[i];

        if (calculated !== expected) {
            console.log("FAILED");
            for (let j = 0; j < test[i][0].length; j++) {
                console.log(test[i][0][j]);
            }
            console.log("Should Return: " + expected);
            console.log("But Returns: " + calculated);
        }

        assert.equal(calculated, expected);
    }
}