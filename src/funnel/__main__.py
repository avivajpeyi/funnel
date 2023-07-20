"""Command-line interface."""
import click


@click.command()
@click.version_option()
def main() -> None:
    """Funnel."""


if __name__ == "__main__":
    main(prog_name="funnel")  # pragma: no cover
