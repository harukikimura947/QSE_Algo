using BenchmarkTools
using LinearAlgebra
using Random
using Profile
using ProfileView
using StatProfilerHTML
using Printf;

println("Before allocation: ", Base.gc_bytes())

# メモリ割り当てを伴う操作
x = rand(1000, 1000);

println("After allocation: ", Base.gc_bytes())

x = nothing

# 明示的にガベージコレクションを実行
GC.gc()

println("After GC: ", Base.gc_bytes())