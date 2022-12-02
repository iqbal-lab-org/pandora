# Warning

The `musl` binary is currently much slower than the `glibc` binary. Do not build this.

# Usage

```
cd <pandora_folder>/scripts/portable_binary_builder/musl
./build_portable_binary.sh
```

The portable binary will be in `<pandora_folder>/pandora-linux-precompiled-musl`

This binary is statically linked to `musl`.
