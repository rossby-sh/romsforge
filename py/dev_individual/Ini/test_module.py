import yaml


class ConfigObject:
    def __init__(self, **entries):
        for key, value in entries.items():
            if isinstance(value, dict):
                value = ConfigObject(**value)
            self.__dict__[key] = value

    def __getitem__(self, key):
        return self.__dict__[key]

    def __repr__(self, indent=0):
        lines = []
        pad = '  ' * indent
        for key, value in self.__dict__.items():
            if isinstance(value, ConfigObject):
                lines.append(f"{pad}{key}:")
                lines.append(value.__repr__(indent + 1))
            else:
                lines.append(f"{pad}{key}: {value}")
        return '\n'.join(lines)

    def to_dict(self):
        out = {}
        for key, value in self.__dict__.items():
            if isinstance(value, ConfigObject):
                out[key] = value.to_dict()
            else:
                out[key] = value
        return out

def parse_config(path="config.yaml") -> ConfigObject:
    with open(path) as f:
        raw = yaml.safe_load(f)
    return ConfigObject(**raw)



cfg = parse_config("./config.yaml")

# dot-access
print(cfg.ininame)
print(cfg.sigma_coord.layer_n)

# dict-access도 가능
print(cfg["sigma_coord"].theta_s)
print(cfg["ogcm_var_name"]["temperature"])
print(cfg.to_dict())
print('!!!')
print(cfg)














