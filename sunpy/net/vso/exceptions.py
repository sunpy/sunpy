class NoData(Exception):
    """
    Risen for callbacks of VSOClient that are unable to supply information for the request.
    """


class DownloadFailed(Exception):
    pass


class MissingInformation(Exception):
    pass


class UnknownMethod(Exception):
    pass


class MultipleChoices(Exception):
    pass


class UnknownVersion(Exception):
    pass


class UnknownStatus(Exception):
    pass
