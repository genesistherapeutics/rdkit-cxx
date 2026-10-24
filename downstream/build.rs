fn main() {
    println!(
        "cargo::rustc-link-arg=-Wl,-rpath,{}",
        std::env::var("DEP_RDKIT_CXX_LIBDIR").unwrap()
    );
}
