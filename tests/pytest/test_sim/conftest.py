import pytest

def pytest_addoption(parser):
    parser.addoption(
        "--mpath", action="store", help="path to mesh files")

@pytest.fixture
def mpath(request):
    return request.config.getoption("--mpath")
