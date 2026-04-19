---
name: authoring-sf-programs-java
description: Use when authoring a Madagascar sf* main program in Java.
---

## When to use

Choose the Java path when researchers bring an existing Java background or need access
to Java scientific libraries (e.g., the Mines JTK) that have no C/Fortran equivalent.
The file naming convention follows the same pattern as other language bindings:
`M<name>.java` (e.g., `MClip.java` produces the program `sfclip`).

**Note:** Java support requires a JDK at configure time and the `JAVA_HOME` (or
`JAVA_SDK`) environment variable set before running SCons. It is rarely used in
practice — there are no `M*.java` programs in the main Madagascar tree. The files
under `api/java/` (`Input.java`, `Output.java`, etc.) are reference implementations
and the `test/Clip.java` example demonstrates usage. Additionally, the Mines JTK
(`MINESJTK` environment variable) is required for the older API layer under
`api/java/old/`.

For the broader context of what an `sf*` program is, how self-documentation works,
and general conventions that apply across all languages, see the shared skill:
[skills/authoring-sf-programs/](../authoring-sf-programs/SKILL.md).

## Skeleton

A minimal Java program that reads a 3-D float dataset, processes it trace-by-trace,
and writes the result. Based on the `test/Clip.java` reference implementation.

```java
import rsf.RSF;
import rsf.Input;
import rsf.Output;

public class MExample {
    // Load the JNI bridge at class-load time
    static {
        System.loadLibrary("jrsf");
    }

    public static void main(String[] args) {
        // 1. Initialize the RSF parameter system with command-line args
        RSF par = new RSF(args);

        // 2. Open standard input and output RSF files
        Input  in  = new Input("in");
        Output out = new Output("out");

        // 3. Read header dimensions (1-based axis indexing)
        int n1 = in.getN(1);
        int n2 = in.getN(2);
        int n3 = in.getN(3);

        // 4. Copy axis metadata to output
        out.setN(1, n1);  out.setDelta(1, in.getDelta(1));
        out.setN(2, n2);  out.setOrigin(1, in.getOrigin(1));
        out.setN(3, n3);

        // 5. Parse a user parameter (with default fallback)
        float threshold = par.getFloat("threshold", 1.0f);

        // 6. Read, process, write one fast-axis slice at a time
        float[] trace = new float[n1];
        for (int i3 = 0; i3 < n3; i3++) {
            for (int i2 = 0; i2 < n2; i2++) {
                in.read(trace);
                // ... your processing here ...
                out.write(trace);
            }
        }

        // 7. Close files
        in.close();
        out.close();
    }
}
```

The `static { System.loadLibrary("jrsf"); }` block is mandatory — it loads the
SWIG-generated JNI shared library (`libjrsf.jnilib` on macOS, `jrsf.so` on Linux)
bridging Java to the underlying `libdrsf` C library.

## API cheat sheet

### RSF — parameter / initialisation (`api/java/RSF.java`)

| Method | Signature | Notes |
|--------|-----------|-------|
| Constructor | `new RSF(String[] args)` | Prepends `"java"` and calls `sf_init`; pass `main`'s `args` directly |
| getBool | `boolean getBool(String key, boolean fallback)` | Returns the parsed bool; `fallback` not yet wired — underlying `sf_getbool` result is returned |
| getInt | `int getInt(String key, int fallback)` | Returns `fallback` when key absent |
| getFloat | `float getFloat(String key, float fallback)` | Returns `fallback` when key absent |
| getString | `String getString(String key, String fallback)` | Returns `fallback` when key absent or null |

### RSFFile — header access (`api/java/RSFFile.java`)

Abstract base for `Input` and `Output`. All index arguments are **1-based**.

| Method | Signature | Notes |
|--------|-----------|-------|
| getN | `int getN(int axis)` | Reads `n<axis>` from header via `sf_histint` |
| setN | `void setN(int axis, int n)` | Writes `n<axis>` via `sf_putint` |
| getDelta | `float getDelta(int axis)` | Reads `d<axis>` via `sf_histfloat` |
| setDelta | `void setDelta(int axis, float d)` | Writes `d<axis>` via `sf_putfloat` |
| getOrigin | `float getOrigin(int axis)` | Reads `o<axis>` via `sf_histfloat` |
| setOrigin | `void setOrigin(int axis, float o)` | Writes `o<axis>` via `sf_putfloat` |
| getLabel | `String getLabel(int axis)` | Reads `label<axis>` via `sf_histstring` |
| setLabel | `void setLabel(int axis, String s)` | Writes `label<axis>` via `sf_putstring` |
| getUnit | `String getUnit(int axis)` | Reads `unit<axis>` via `sf_histstring` |
| setUnit | `void setUnit(int axis, String s)` | Writes `unit<axis>` via `sf_putstring` |
| close | `void close()` | Calls `sf_fileclose`; always call for both `Input` and `Output` |
| MAX_DIMS | `static final int` | Equals `SF_MAX_DIM` from the generated `m8rConstants` class |

### Input — reading data (`api/java/Input.java`)

Extends `RSFFile`. Delegates to `sf_floatread` via the SWIG-generated `m8r` class.

| Method | Notes |
|--------|-------|
| `new Input("in")` | `"in"` is the standard stdin token; any RSF filename also works |
| `void read(float[] a)` | Read `a.length` floats from the data stream |
| `void read(float[][] a)` | Iterates over outer dimension, calls `read(float[])` |
| `void read(float[][]...a)` | Overloads up to 9-D arrays via recursive delegation |

### Output — writing data (`api/java/Output.java`)

Extends `RSFFile`. Delegates to `sf_floatwrite`. Wraps all writes in
`try/catch (Exception e)`, printing to `System.err` on failure.

| Method | Notes |
|--------|-------|
| `new Output("out")` | `"out"` is the standard stdout token |
| `void write(float[] d)` | Write `d.length` floats to the data stream |
| `void write(float[][] d)` | Iterates and delegates to `write(float[])` |
| `void write(float[][]...d)` | Overloads up to 9-D arrays |

### Error handling

`Output` and `RSFFile` accessors catch `Exception`, print to `System.err`, and return
sensible defaults (0 / 0.0f / `""`). No Java-level `sf_error`; fatal errors propagate from the C layer.

## Build integration

The `api/java/SConstruct` file documents the full build process:

1. **JDK detection** — reads `JAVA_HOME` or `JAVA_SDK` from the SCons environment.
   If neither is set the build exits immediately.
2. **SWIG step** — `m8r.i` is processed by SWIG with `-java -package rsf` to
   produce `m8r_wrap.c` and the generated Java stubs (`m8r.java`, `m8rJNI.java`,
   `m8rConstants.java`, `SWIGTYPE_p_sf_File.java`, etc.).
3. **JNI shared library** — the wrap C file is compiled and linked against
   `libdrsf` to produce `libjrsf.jnilib` (macOS) or `jrsf.so` (Linux).
4. **Java compilation** — `javac` compiles `RSF.java`, `RSFFile.java`,
   `Input.java`, `Output.java`, and the generated stubs into `rsf/*.class`.
5. **JAR packaging** — `rsf.class` files are archived into `rsf.jar` and
   installed to `$RSFROOT/lib/`.

For a user program there is no `UserSconsTargets` Java attribute in the SCons
framework. User programs must replicate the compile-and-jar pattern manually or
use the `project.Java(...)` helper demonstrated in `api/java/test/SConstruct`:

```python
from rsf.proj import *
project.Java('.', 'MExample.java')
Flow('result', 'input MExample.class',
     '%s MExample threshold=0.5' % WhereIs('java'))
End()
```

The JVM must be able to locate both `rsf.jar` and the JNI library at runtime.
Set `CLASSPATH` and `java.library.path` (or `LD_LIBRARY_PATH` / `DYLD_LIBRARY_PATH`)
to point at `$RSFROOT/lib/`.

## Pointers

All files under `api/java/`:

| File | Description |
|------|-------------|
| `README` | Short note on the Mines JTK dependency and `MINESJTK` environment variable |
| `RSF.java` | Top-level class: calls `sf_init` and exposes `getInt/getFloat/getString/getBool` for command-line parameter parsing |
| `RSFFile.java` | Abstract base class for RSF files: header getters/setters (`getN`, `getDelta`, `getOrigin`, `getLabel`, `getUnit`) and `close()` |
| `Input.java` | Concrete read-only RSF file: `read(float[])` overloaded up to 9 dimensions |
| `Output.java` | Concrete write-only RSF file: `write(float[])` overloaded up to 9 dimensions |
| `m8r.i` | SWIG interface file that generates the JNI bridge between Java and the Madagascar C API (`libdrsf`) |
| `SConstruct` | SCons build script: JDK detection, SWIG invocation, `javac` compilation, JAR creation, and install to `$RSFROOT/lib/` |
| `test/Clip.java` | Reference implementation of a clipping program; shows the complete pattern: `RSF`, `Input`, `Output`, header reads, trace loop, `read`/`write`, `close` |
| `test/SConstruct` | Minimal SCons flow using `project.Java()` and `WhereIs('java')` to invoke the compiled class |
| `old/Header.java` | Legacy header class (Mines JTK-based); superseded by `RSFFile.java` |
| `old/Par.java` | Legacy parameter-parsing class; superseded by `RSF.java` |
| `old/Reader.java` | Legacy reader; superseded by `Input.java` |
| `old/Writer.java` | Legacy writer; superseded by `Output.java` |

## Shared conventions

Self-documentation, program naming (`sf*` / `M*.java`), data-type conventions,
axis ordering (n1 is the fast axis), and the overall SCons flow are shared across
all Madagascar language bindings. See:

[skills/authoring-sf-programs/](../authoring-sf-programs/SKILL.md)
