# PalabosDev

## install the dependency packages by vcpkg

`vcpkg install --feature-flags=binarycaching,manifests`

## Generate core file when meet a core dump

`ulimit -c`

if result is '0', no coredump file will be generated

`ulimit -c unlimited`

echo it into .zshrc to make it permanent.

`ulimit -c unlimited`

## Using gdb to debug core dump

`gdb <your program> <core file>`
