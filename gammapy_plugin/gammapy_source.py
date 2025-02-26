import collections
import numpy as np
from astromodels.sources.source import Source, SourceType
from astromodels.core.tree import Node
from astromodels.core.parameter import Parameter
from astromodels.utils.pretty_list import dict_to_list
from gammapy.modeling.models.core import ModelBase

# TODO docs
# TODO typ hinting


class GammapySource(Source, Node):
    """
    A dummy source for using a GammaPy Model
    Only works in the case of it being directly assigned to

    """

    def __init__(self, name, gammapy_model: ModelBase):
        self._gammapy_model = gammapy_model
        Source.__init__(self, [], SourceType.EXTENDED_SOURCE)
        Node.__init__(self, name)
        self._get_parameters()
        spectrum_node = Node("spectrum")
        self._add_child(spectrum_node)

    def _get_parameters(self):
        self._parameters = collections.OrderedDict()
        tp = self._gammapy_model.parameters.to_dict()
        parameters = collections.OrderedDict()
        for i in range(len(tp)):
            ttp = tp[i]
            norm = False
            if ttp["name"] == "norm":
                # TODO check if this covers all cases
                norm = True
            min_v = ttp["min"]
            max_v = ttp["max"]
            if np.isnan(min_v):
                min_v = None
            if np.isnan(max_v):
                max_v = None
            parameters[ttp["name"]] = Parameter(
                name=ttp["name"],
                value=ttp["value"],
                min_value=min_v,
                max_value=max_v,
                is_normalization=norm,
                free=not ttp["frozen"],
                unit=ttp["unit"],
                desc=f"Gammapy Model Parameter {ttp['name']}",
            )
        for child_name, child in parameters.items():
            self._parameters[child_name] = child
            self._add_child(child)

    @property
    def parameters(self):
        return self._parameters

    @property
    def free_parameters(self):
        return [p for k, p in self._parameters.items() if p.free]

    @property
    def has_free_parameters(self):
        ret = False
        for k, v in self._parameters.items():
            if v.free:
                ret = True
                break
        return ret

    @property
    def gammapy_model(self):
        return self._gammapy_model

    @property
    def linked_parameters(self):
        return self._wrapped.linked_parameters

    def _repr__base(self, rich_output=False):

        repr_dict = collections.OrderedDict()
        key = "%s (gammapy source)" % self.name
        repr_dict[key] = collections.OrderedDict()
        for component_name, component in list(self.components.items()):
            repr_dict[key][component_name] = component.to_dict(minimal=True)

        return dict_to_list(repr_dict, rich_output)
