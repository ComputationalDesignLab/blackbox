import inspect
from importlib import import_module
from typing import get_args, get_origin, Literal

from docutils import nodes
from docutils.statemachine import StringList
from sphinx.util.docutils import SphinxDirective

from pydantic import BaseModel

class PydanticOptions(SphinxDirective):
    """
        Custom Sphinx directive for reading a pydantic class and
        converting it to appropriate rendering
    """

    required_arguments = 1
    optional_arguments = 0
    has_content = False

    def run(self):

        module_name, class_name = self.arguments[0].rsplit(".", 1)

        module = import_module(module_name)
        model = getattr(module, class_name)

        if not issubclass(model, BaseModel):
            raise TypeError(
                f"{self.arguments[0]} must inherit from BaseModel"
            )

        # --------------------------------------------------
        # Generate RST in memory
        # --------------------------------------------------

        lines = []

        for name, field in model.model_fields.items():

            lines.extend(
                self.make_option_rst(
                    name,
                    field,
                )
            )

            # Blank line between data directives
            lines.append("")

        # --------------------------------------------------
        # Let Sphinx parse the generated RST
        # --------------------------------------------------

        container = nodes.container(
            classes=["pydantic-options"]
        )

        rst = StringList(
            lines,
            source=self.state.document.current_source,
        )

        self.state.nested_parse(
            rst,
            self.content_offset,
            container,
        )

        return [container]

    # def run(self):

    #     module_name, class_name = self.arguments[0].rsplit(".", 1)

    #     module = import_module(module_name)
    #     model = getattr(module, class_name)

    #     if not issubclass(model, BaseModel):
    #         raise TypeError(
    #             f"{self.arguments[0]} must inherit from BaseModel"
    #         )

    #     fields = model.model_fields

    #     groups = {}

    #     for name, field in fields.items():

    #         extra = field.json_schema_extra or {}

    #         group = extra.get(
    #             "group",
    #             "Other arguments",
    #         )

    #         groups.setdefault(group, []).append(
    #             (name, field)
    #         )

    #     container = nodes.container(
    #         classes=["pydantic-options"]
    #     )

    #     for group_name, group_fields in groups.items():

    #         # Group heading
    #         container += nodes.rubric(
    #             "",
    #             group_name,
    #             classes=["pydantic-options-group"],
    #         )

    #         # Options
    #         lines = []

    #         for name, field in group_fields:

    #             lines.extend(
    #                 self.make_option_rst(
    #                     name,
    #                     field,
    #                 )
    #             )

    #             lines.append("")

    #         rst = StringList(``
    #             lines,
    #             source=self.state.document.current_source,
    #         )

    #         self.state.nested_parse(
    #             rst,
    #             self.content_offset,
    #             container,
    #         )

    #     return [container]

    def make_option_rst(self, name, field):

        lines = []

        # --------------------------------------------------
        # Data directive
        # --------------------------------------------------

        lines.append(
            f".. py:data:: {name}"
        )

        # Type
        lines.append(
            f"   :type: {self.format_type(field.annotation)}"
        )

        # Default
        if not field.is_required():
            lines.append(
                f"   :value: {self.format_default(field)}"
            )

        # --------------------------------------------------
        # Description
        # --------------------------------------------------

        if field.description:

            lines.append("")

            description = inspect.cleandoc(
                field.description
            )

            for line in description.splitlines():

                if line:
                    lines.append(
                        f"   {line}"
                    )
                else:
                    lines.append("")

        return lines

    @staticmethod
    def format_type(annotation):

        origin = get_origin(annotation)
        args = get_args(annotation)

        # ----------------------------------------
        # list[...]
        # ----------------------------------------

        if origin is list:

            # For documentation, display simply as "list"
            return "list"

        # ----------------------------------------
        # Literal[...]
        # ----------------------------------------

        if origin is Literal:

            if not args:
                return "Literal"

            # Infer the underlying type from the
            # literal values.
            first_value = args[0]

            if isinstance(first_value, bool):
                return "bool"

            if isinstance(first_value, int):
                return "int"

            if isinstance(first_value, float):
                return "float"

            if isinstance(first_value, str):
                return "str"

            return type(first_value).__name__

        # ----------------------------------------
        # dict[...]
        # ----------------------------------------

        if origin is dict:

            if len(args) == 2:
                return (
                    "dict["
                    f"{PydanticOptions.format_type(args[0])}, "
                    f"{PydanticOptions.format_type(args[1])}"
                    "]"
                )

            return "dict"

        # ----------------------------------------
        # Normal types
        # ----------------------------------------

        if hasattr(annotation, "__name__"):
            return annotation.__name__

        return str(annotation).replace(
            "typing.",
            "",
        )

    @staticmethod
    def format_default(field):

        if field.default_factory is not None:
            default = field.default_factory()
        else:
            default = field.default

        if isinstance(default, str):
            return f"'{default}'"

        return repr(default)

def setup(app):

    app.add_directive(
        "pydantic-options",
        PydanticOptions,
    )

    return {
        "version": "0.1",
        "parallel_read_safe": True,
    }
