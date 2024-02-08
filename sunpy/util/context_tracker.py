__all__ = ["ContextTracker"]

class ContextTracker:
    """
    A class to keep track of the current context of the code.
    """
    def __init__(self):
        self._active_contexts = []

    def enter_context(self, context):
        self._active_contexts.append(context)

    def exit_context(self, context):
        self._active_contexts.remove(context)

    def is_active(self,context):
        return context in self._active_contexts


global_context_tracker = ContextTracker()
