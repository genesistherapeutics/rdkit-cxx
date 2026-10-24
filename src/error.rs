/// Unwraps the result of a bridge call that RDKit is not expected to fail.
///
/// An exception here is a bug in the bindings, not a user error (E1).
pub(crate) fn expect_ok<T>(result: Result<T, cxx::Exception>) -> T {
    match result {
        Ok(value) => value,
        Err(exception) => panic!("unexpected RDKit exception: {}", exception.what()),
    }
}
