use anyhow::{Context, Result};
use primitive_types::U256;
use starkom_bluesky::Scalar;

/// Converts a BlueSky scalar to U256.
pub fn scalar_to_u256(value: Scalar) -> U256 {
    U256::from_little_endian(&value.to_little_endian())
}

/// Converts a U256 value to a BlueSky scalar, failing if the value is outside the BlueSky range.
pub fn u256_to_scalar(value: U256) -> Result<Scalar> {
    Scalar::from_repr_vartime(&value.to_little_endian()).context("invalid BlueSky scalar")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse_scalar(s: &'static str) -> Scalar {
        s.parse().unwrap()
    }

    #[test]
    fn test_scalar_to_u256() {
        assert_eq!(
            scalar_to_u256(parse_scalar(
                "0x9ff20c13ccb8a61ced7558c8e10964efd5ee3557d3a2bc0dfb83662950fc85f"
            )),
            "0x9ff20c13ccb8a61ced7558c8e10964efd5ee3557d3a2bc0dfb83662950fc85f"
                .parse()
                .unwrap()
        );
    }

    #[test]
    fn test_u256_to_scalar() {
        assert_eq!(
            u256_to_scalar(
                "0x18d82aec545e64ec800bfd5d81baed36fa8c3ea2fdf5514256eb5bf312613a8e"
                    .parse()
                    .unwrap()
            )
            .unwrap(),
            parse_scalar("0x18d82aec545e64ec800bfd5d81baed36fa8c3ea2fdf5514256eb5bf312613a8e")
        );
    }
}
