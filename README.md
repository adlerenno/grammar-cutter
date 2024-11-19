# grammar-cutter

**Autor:** EA, DK

grammar-cutter is a C-library using [SERP]() to generate RePair grammars and to extract substrings as compressed grammars from a compressed grammar.

The library should run on Linux and macOS, it is not supported on Windows.

If you have questions regarding the implementation, feel free to contact adlerenno.

## Dependencies

- none

## Build

To build, you need `cmake` and you can compile the library by:

```bash
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DOPTIMIZE_FOR_NATIVE=on ..
make
```

The following parameters can be passed to CMake:

- `-DCMAKE_BUILD_TYPE=Release` activates compiler optimizations
- `-DOPTIMIZE_FOR_NATIVE=on` activates optimized processor functions, for example the more efficient `popcnt`-variants
- `-DCLI=on` activates the compilation of the command-line-tool.
- `-DWEB_SERVICE=on` activates the compilation of the web-service. Needs the command-line-tool.

The library will be in the build-directory as "libgrammar-cutter.1.0.0.dylib" (macOS) or "libgrammar-cutter.so.1.0.0" (Linux).
The command-line-tool is in the build-directory as well and is called "libgrammar-cutter-cli".

## Command-Line-Tool

With the command-line-tool "libgrammar-cutter-cli" you can compress graphs and search within them.
The following help can also be obtained by `./libgrammar-cutter-cli --help`.

```
Usage: 
	./grammar-cutter-cli [-i <input_file>] [-o <output_file>] [-f <from_position>] [-t <to_position>] 
or
	./grammar-cutter-cli -h
```

You can use ```-i``` and ```-o``` together to compress a file or ```-i```, ```-f``` and ```-t``` to query an interval.

## Library

TODO
