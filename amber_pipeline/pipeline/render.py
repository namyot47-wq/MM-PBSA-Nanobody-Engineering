#pipeline/render.py 

from jinja2 import Environment, FileSystemLoader
def render_input(template_name: str, context: dict, output_path: str):
    env = Environment(loader=FileSystemLoader("equilibrium_templates"))
    Path(out_path).write_text(template.render(**context))
    