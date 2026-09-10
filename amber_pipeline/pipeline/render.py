#pipeline/render.py
from pathlib import Path
from jinja2 import Environment, FileSystemLoader

def render_input(template_name: str, context: dict, output_path: str):
    env = Environment(loader=FileSystemLoader("equilibrium_templates"))
    template = env.get_template(template_name)
    Path(output_path).write_text(template.render(**context))
