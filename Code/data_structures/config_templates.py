from pathlib import Path
from typing import Optional, Union

from pydantic import BaseModel


class GisFile(BaseModel):
    column_mapping: dict
    column_selection: Optional[dict]
    name: str
    path: Union[str, list, dict, Path]


class FMGisFile(GisFile):
    pass


class RRGisFile(GisFile):
    pass
