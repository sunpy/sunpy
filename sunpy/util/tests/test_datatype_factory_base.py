from __future__ import absolute_import, division, print_function

import pytest

from sunpy.util.datatype_factory_base import BasicRegistrationFactory
from sunpy.util.datatype_factory_base import NoMatchError
from sunpy.util.datatype_factory_base import MultipleMatchError
from sunpy.util.datatype_factory_base import ValidationFunctionError


class BaseWidget(object):
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


class TestBasicRegistrationFactory(object):

    def test_default_factory(self):

        DefaultFactory = BasicRegistrationFactory()

        DefaultFactory.register(DefaultWidget, is_default=True)
        assert DefaultFactory.default_widget_type == DefaultWidget

        DefaultFactory.register(StandardWidget)
        DefaultFactory.register(FancyWidget)
        DefaultFactory.register(ExternallyValidatedWidget,
                                validation_function=external_validation_function)

        assert type(DefaultFactory()) is DefaultWidget
        assert type(DefaultFactory(style='standard')) is StandardWidget
        assert type(DefaultFactory(style='fancy')) is DefaultWidget
        assert type(DefaultFactory(style='fancy', feature="present")) is FancyWidget
        assert type(DefaultFactory(style='external')) is ExternallyValidatedWidget

        with pytest.raises(ValidationFunctionError):
            DefaultFactory.register(UnvalidatedWidget)

        with pytest.raises(ValidationFunctionError):
            DefaultFactory.register(MissingClassMethodWidget)

        DefaultFactory.unregister(StandardWidget)
        assert type(DefaultFactory(style='standard')) is not StandardWidget

    def test_validation_fun_not_callable(self):
        TestFactory = BasicRegistrationFactory()

        with pytest.raises(AttributeError):
            TestFactory.register(StandardWidget, validation_function='not_callable')

    def test_no_default_factory(self):

        NoDefaultFactory = BasicRegistrationFactory()

        NoDefaultFactory.register(StandardWidget)
        NoDefaultFactory.register(FancyWidget)

        with pytest.raises(NoMatchError):
            NoDefaultFactory()

        # Raises because all requirements are not met for FancyWidget and no
        # default is present.
        with pytest.raises(NoMatchError):
            NoDefaultFactory(style='fancy')

        assert type(NoDefaultFactory(style='standard')) is StandardWidget
        assert type(NoDefaultFactory(style='fancy', feature='present')) is FancyWidget

    def test_with_external_registry(self):
        external_registry = {}

        FactoryWithExternalRegistry = \
            BasicRegistrationFactory(registry=external_registry)

        assert len(external_registry) == 0

        FactoryWithExternalRegistry.register(StandardWidget)
        assert type(FactoryWithExternalRegistry(style='standard')) is StandardWidget

        # Ensure the 'external_registry' is being populated see #1988
        assert len(external_registry) == 1

    def test_multiple_match_factory(self):

        MultipleMatchFactory = BasicRegistrationFactory()

        MultipleMatchFactory.register(StandardWidget)
        MultipleMatchFactory.register(DuplicateStandardWidget)

        with pytest.raises(MultipleMatchError):
            MultipleMatchFactory(style='standard')

    def test_extra_validation_factory(self):
        ExtraValidationFactory = \
            BasicRegistrationFactory(
                additional_validation_functions=['different_validation_function'])

        ExtraValidationFactory.register(DifferentValidationWidget)

        assert type(ExtraValidationFactory(style='different')) is DifferentValidationWidget

        with pytest.raises(ValidationFunctionError):
            ExtraValidationFactory.register(MissingClassMethodDifferentValidationWidget)
