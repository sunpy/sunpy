import pytest

from sunpy.util.datatype_factory_base import (
    BasicRegistrationFactory,
    MultipleMatchError,
    NoMatchError,
    ValidationFunctionError,
)


class BaseWidget:
    def __init__(self, *args, **kwargs):
        pass


class DefaultWidget(BaseWidget):
    pass


class StandardWidget(BaseWidget):
    @classmethod
    def _factory_validation_function(cls, *args, **kwargs):
        return kwargs.get('style') == 'standard'


class DuplicateStandardWidget(BaseWidget):
    @classmethod
    def _factory_validation_function(cls, *args, **kwargs):
        return kwargs.get('style') == 'standard'


class FancyWidget(BaseWidget):
    @classmethod
    def _factory_validation_function(cls, *args, **kwargs):
        return kwargs.get('style') == 'fancy' and 'feature' in kwargs


class ExternallyValidatedWidget(BaseWidget):
    pass


def external_validation_function(*args, **kwargs):
    return kwargs.get('style') == 'external'


class UnvalidatedWidget(BaseWidget):
    pass


class MissingClassMethodWidget(BaseWidget):
    def _factory_validation_function(cls, *args, **kwargs):
        return kwargs.get('style') == 'missing'


class DifferentValidationWidget(BaseWidget):
    @classmethod
    def different_validation_function(cls, *args, **kwargs):
        return kwargs.get('style') == 'different'


class MissingClassMethodDifferentValidationWidget(BaseWidget):
    def different_validation_function(cls, *args, **kwargs):
        return kwargs.get('style') == 'missing-different'


def test_default_factory():
    DefaultFactory = BasicRegistrationFactory()

    DefaultFactory.register(DefaultWidget, is_default=True)
    assert DefaultFactory.default_widget_type == DefaultWidget

    DefaultFactory.register(StandardWidget)
    DefaultFactory.register(FancyWidget)
    DefaultFactory.register(ExternallyValidatedWidget,
                            validation_function=external_validation_function)

    assert isinstance(DefaultFactory(), DefaultWidget)
    assert isinstance(DefaultFactory(style='standard'), StandardWidget)
    assert isinstance(DefaultFactory(style='fancy'), DefaultWidget)
    assert isinstance(DefaultFactory(style='fancy', feature="present"), FancyWidget)
    assert isinstance(DefaultFactory(style='external'), ExternallyValidatedWidget)

    with pytest.raises(ValidationFunctionError):
        DefaultFactory.register(UnvalidatedWidget)

    with pytest.raises(ValidationFunctionError):
        DefaultFactory.register(MissingClassMethodWidget)

    DefaultFactory.unregister(StandardWidget)
    assert not isinstance(DefaultFactory(style='standard'), StandardWidget)


def test_validation_fun_not_callable():
    TestFactory = BasicRegistrationFactory()

    with pytest.raises(AttributeError):
        TestFactory.register(StandardWidget, validation_function='not_callable')


def test_no_default_factory():
    NoDefaultFactory = BasicRegistrationFactory()
    NoDefaultFactory.register(StandardWidget)
    NoDefaultFactory.register(FancyWidget)

    with pytest.raises(NoMatchError):
        NoDefaultFactory()

    # Raises because all requirements are not met for FancyWidget and no
    # default is present.
    with pytest.raises(NoMatchError):
        NoDefaultFactory(style='fancy')

    assert isinstance(NoDefaultFactory(style='standard'), StandardWidget)
    assert isinstance(NoDefaultFactory(style='fancy', feature='present'), FancyWidget)


def test_with_external_registry():
    external_registry = {}

    FactoryWithExternalRegistry = BasicRegistrationFactory(registry=external_registry)
    assert len(external_registry) == 0

    FactoryWithExternalRegistry.register(StandardWidget)
    assert isinstance(FactoryWithExternalRegistry(style='standard'), StandardWidget)

    # Ensure the 'external_registry' is being populated see #1988
    assert len(external_registry) == 1


def test_multiple_match_factory():
    MultipleMatchFactory = BasicRegistrationFactory()
    MultipleMatchFactory.register(StandardWidget)
    MultipleMatchFactory.register(DuplicateStandardWidget)

    with pytest.raises(MultipleMatchError):
        MultipleMatchFactory(style='standard')


def test_extra_validation_factory():
    ExtraValidationFactory = BasicRegistrationFactory(
        additional_validation_functions=['different_validation_function'])
    ExtraValidationFactory.register(DifferentValidationWidget)
    assert isinstance(ExtraValidationFactory(style='different'), DifferentValidationWidget)

    with pytest.raises(ValidationFunctionError):
        ExtraValidationFactory.register(MissingClassMethodDifferentValidationWidget)
