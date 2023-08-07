import click
from pathlib import Path
from miqssr.model_build import build_model


main = click.Group()


@main.command(name='build_model')
@click.option(
    "--config",
    "config",
    required=True,
    help="Path to the config YAML",
    type=click.Path(exists=True, path_type=Path),
)
def build_model_cli(config):
    build_model(config)


if __name__ == '__main__':
    main()

# python GSLRetro/interfaces/cli.py download_data



