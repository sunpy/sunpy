from sunpy.util.datatype_factory_base import BasicRegistrationFactory
from sunpy.net.attr import *
from sunpy.net.vso.attrs import *
from sunpy.util.datatype_factory_base import NoMatchError
from sunpy.util.datatype_factory_base import MultipleMatchError

__all__=['UnifiedDownloader']

qwalker = AttrWalker()

@qwalker.add_creator(AttrAnd)
def _create(wlk,query,dobj):
    
    qresponseobj,qclient = dobj._get_registered_widget(*query.attrs)
    return [(qresponseobj,qclient)]


@qwalker.add_creator(AttrOr)
def _create(wlk,query,dobj):    
    
    qblocks = []
    for iattr in query.attrs:
        qblocks.extend(wlk.create(iattr,dobj))
    
    return qblocks


class UnifiedDownloaderFactory(BasicRegistrationFactory):


    def query(self,*query):
        
        query = and_(*query)
        return qwalker.create(query,self)
    
    def __call__(self,*args,**kwargs):
        pass


    def _check_registered_widgets(self,*args,**kwargs):
    
        candidate_widget_types = list()
	for key in self.registry:
	    
	    if self.registry[key](*args):
	        candidate_widget_types.append(key)
            
	n_matches = len(candidate_widget_types)
	if n_matches == 0:
	    if self.default_widget_type is None:
	        raise NoMatchError("Query {0} can not be handled in its current form".format(args))
	    else:
	        return  [self.default_widget_type]
	elif n_matches > 1:
	    raise MultipleMatchError("Too many candidates clients can service your query {0}".format(args))

        return candidate_widget_types

    def _get_registered_widget(self,*args,**kwargs):

        candidate_widget_types = self._check_registered_widgets(*args)
	tmpclient = candidate_widget_types[0]()
	return tmpclient.query(*args),tmpclient


UnifiedDownloader = UnifiedDownloaderFactory(additional_validation_functions = ['_can_handle_query'])
