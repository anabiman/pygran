import pytest


def pytest_addoption(parser):
    parser.addoption("--trajf", action="store", help="trajectory filename")


@pytest.fixture
def trajf(request):
    return request.config.getoption("--trajf")
