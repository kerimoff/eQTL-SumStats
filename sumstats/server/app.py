import numpy as np
import simplejson
import json
from flask import Flask, url_for, request, abort, make_response

import argparse
from collections import OrderedDict
from urllib.parse import unquote
import sumstats.explorer as ex
import sumstats.search as search
from config import properties
from sumstats.common_constants import *
from sumstats.utils.interval import *
from sumstats.server.error_classes import *
from sumstats.errors.error_classes import *

app = Flask(__name__)
app.url_map.strict_slashes = False


@app.errorhandler(APIException)
def handle_custom_api_exception(error):
    response = simplejson.dumps(error.to_dict())
    return make_response(response, error.status_code)


@app.errorhandler(404)
def not_found(error):
    return make_response(simplejson.dumps({'error': 'Page not found'}), 404)


def _set_properties():
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-config', help='The configuration file to use instead of default')
    args = argparser.parse_args()
    if args.config is not None:
        with open(args.config) as config:
            props = json.load(config)
            properties.output_path = props["output_path"]
            properties.gwas_study_location = props["gwas_study_location"]
            properties.input_path = props["input_path"]
            properties.ols_terms_location = props["ols_terms_location"]


def _reconstruct_pvalue(mantissa_dset, exp_dset):
    pval_array = np.empty(len(mantissa_dset), dtype=vlen_dtype)
    for index, mantissa in enumerate(mantissa_dset):
        pval_array[index] = (str(mantissa) + "e" + str(exp_dset[index]))
    return pval_array.tolist()


def _create_href(method_name, params=None):
    params = params or {}
    return {'href': unquote(
        url_for(method_name, **params, _external=True)
    )}


def _create_next_links(method_name, start, size, index_marker, size_retrieved, params=None):
    params = params or {}
    prev = max(0, start - size)
    start_new = start + index_marker

    response = {'self': _create_href(method_name=method_name, params=params)}
    params['start'] = 0
    params['size'] = size
    response['first'] = _create_href(method_name=method_name, params=params)
    params['start'] = prev
    response['prev'] = _create_href(method_name=method_name, params=params)

    if size_retrieved == size:
        params['start'] = start_new
        response['next'] = _create_href(method_name=method_name, params=params)

    return response


def _create_associations_response(method_name, start, size, index_marker, data_dict, params):
    return {'_embedded': {'associations': data_dict}, '_links': _create_next_links(
        method_name=method_name, start=start, size=size, index_marker=index_marker,
        size_retrieved=len(data_dict),
        params=params
    )}


def _create_study_info_for_trait(studies, trait):
    study_list = []
    for study in studies:
        study_info = {'study': study, 'trait': trait,
                      '_links': {'self': _create_href(method_name='get_trait_study_assocs',
                                                      params={'trait': trait, 'study': study}),
                                 'trait': _create_href(method_name='get_trait_assocs', params={'trait': trait})}}

        study_info['_links']['gwas_catalog'] = {'href': str(properties.gwas_study_location + study)}
        study_info['_links']['ols'] = {'href': str(properties.ols_terms_location + trait)}
        study_list.append(study_info)
    return study_list


def _get_array_to_display(datasets, variant=None):
    if datasets is None: return {}
    if len(datasets[REFERENCE_DSET]) <= 0: return {}

    mantissa_dset = datasets.pop(MANTISSA_DSET)
    exponent_dset = datasets.pop(EXP_DSET)
    datasets[PVAL_DSET] = _reconstruct_pvalue(mantissa_dset=mantissa_dset, exp_dset=exponent_dset)

    data_dict = {}
    length = len(datasets[PVAL_DSET])
    for index in range(length):
        # elements are numpy types, they need to be python types to be json serializable
        element_info = {dset: np.asscalar(np.array(dataset[index])) for dset, dataset in datasets.items()}

        if variant is not None:
            var = variant
        else:
            var = datasets[SNP_DSET][index]

        element_info['_links'] = {'self': _create_href(method_name='get_variants',
                                                   params={'variant': var, 'study_accession': datasets[STUDY_DSET][index]})}
        element_info['_links']['variant'] = _create_href(method_name='get_variants',
                                                         params={'variant': var})
        element_info['_links']['study'] = _create_href(method_name='get_studies',
                                                       params={'looking_for': datasets[STUDY_DSET][index]})

        data_dict[index] = element_info
    return data_dict


def _retrieve_endpoint_arguments(args, argument_name, value_if_empty=None):
    try:
        argument = args[argument_name]
    except KeyError:
        argument = value_if_empty
    return argument


def _get_interval(value, interval):
    if value is not None:
        if ":" not in value:
            value = value + ":" + value
        return interval().set_string_tuple(value)
    return None


def _get_basic_arguments(args):
    start = int(_retrieve_endpoint_arguments(args, "start", 0))
    size = int(_retrieve_endpoint_arguments(args, "size", 20))
    pval = _retrieve_endpoint_arguments(args, "p-value")
    pval_interval = _get_interval(pval, FloatInterval)
    return start, size, pval, pval_interval


@app.route('/')
def root():
    # args_chromosome = {'p-value': '{lower:upper}', 'bp': '{lower:upper}', 'study_accession': '{study_accession}'}
    # args_variant = {'p-value': '{lower:upper}', 'study_accession': '{study_accession}'}
    response = {
        '_links': {
            'associations': _create_href(method_name='get_assocs', params={'p-value': '{lower:upper}'}),
            'traits': _create_href(method_name='get_traits'),
            'studies': _create_href(method_name='get_studies'),
            'chromosomes': _create_href(method_name='get_chromosomes'),
            'variant': _create_href(method_name='get_variants', params={'variant': '{variant_id}'})
        }
    }
    return simplejson.dumps(OrderedDict(response))


@app.route('/associations')
def get_assocs():
    args = request.args.to_dict()
    start, size, pval, pval_interval = _get_basic_arguments(args)

    searcher = search.Search(properties.output_path)

    # try:
    datasets, index_marker = searcher.search_all_assocs(start=start, size=size, pval_interval=pval_interval)

    data_dict = _get_array_to_display(datasets)
    params = {'p-value': pval}
    response = _create_associations_response(method_name='get_assocs', start=start, size=size, index_marker=index_marker,
                                             data_dict=data_dict, params=params)

    return simplejson.dumps(OrderedDict(response), ignore_nan=True)

    # except ValueError:
    #     abort(404)


@app.route('/traits')
def get_traits():
    explorer = ex.Explorer(properties.output_path)
    trait_list = []
    traits = explorer.get_list_of_traits()
    for trait in traits:
        trait_info = {'trait': trait,
                      '_links': {'self': _create_href(method_name='get_trait_assocs', params={'trait': trait})}}
        trait_info['_links']['studies'] = _create_href(method_name='get_studies_for_trait', params={'trait': trait})
        trait_info['_links']['ols'] = {'href': str(properties.ols_terms_location + trait)}

        trait_list.append(trait_info)

    response = {'_embedded': {'traits': trait_list}}
    return simplejson.dumps(OrderedDict(response))


@app.route('/traits/<string:trait>')
def get_trait_assocs(trait):
    args = request.args.to_dict()
    start, size, pval, pval_interval = _get_basic_arguments(args)

    searcher = search.Search(properties.output_path)

    try:
        datasets, index_marker = searcher.search_trait(trait=trait, start=start, size=size, pval_interval=pval_interval)

        data_dict = _get_array_to_display(datasets)
        params = {'trait': trait, 'p-value': pval}
        response = _create_associations_response(method_name='get_trait_assocs', start=start, size=size, index_marker=index_marker,
                                                 data_dict=data_dict, params=params)
        response['_links']['studies'] = _create_href(method_name='get_studies_for_trait', params={'trait': trait})

        return simplejson.dumps(OrderedDict(response), ignore_nan=True)

    except NotFoundError as error:
        raise RequestedNotFound(str(error))


@app.route('/studies')
@app.route('/studies/<looking_for>')
def get_studies(looking_for=None):
    explorer = ex.Explorer(properties.output_path)
    study_list = []
    trait_studies = []
    if looking_for is not None:
        try:
            trait_studies.append(explorer.get_info_on_study(looking_for))
        except NotFoundError as error:
            raise RequestedNotFound(str(error))
    else:
        trait_studies = explorer.get_list_of_studies()
    for trait_study in trait_studies:
        trait = trait_study.split(":")[0]
        study = trait_study.split(":")[1]

        study_list.append(_create_study_info_for_trait([study], trait))

    response = {'_embedded': {'studies': study_list}}
    return simplejson.dumps(OrderedDict(response))


@app.route('/traits/<string:trait>/studies')
def get_studies_for_trait(trait):
    explorer = ex.Explorer(properties.output_path)
    try:
        studies = explorer.get_list_of_studies_for_trait(trait)
        study_list = _create_study_info_for_trait(studies, trait)
        response = {'_embedded': {'studies': study_list}}
        return simplejson.dumps(OrderedDict(response))
    except NotFoundError as error:
        raise RequestedNotFound(str(error))


@app.route('/traits/<string:trait>/studies/<string:study>')
def get_trait_study_assocs(trait, study):
    args = request.args.to_dict()
    start, size, pval, pval_interval = _get_basic_arguments(args)

    searcher = search.Search(properties.output_path)

    try:
        datasets, index_marker = searcher.search_study(trait=trait, study=study,
                                                       start=start, size=size, pval_interval=pval_interval)

        data_dict = _get_array_to_display(datasets)
        params = {'trait': trait, 'study': study, 'p-value': pval}
        response = _create_associations_response(method_name='get_trait_study_assocs', start=start, size=size, index_marker=index_marker,
                                                 data_dict=data_dict, params=params)

        return simplejson.dumps(OrderedDict(response), ignore_nan=True)

    except NotFoundError as error:
        raise RequestedNotFound(str(error))


@app.route('/chromosomes')
def get_chromosomes():
    chromosomes_list = []
    for chromosome in range(1, 24):
        chromosome_info = {'chromosome': chromosome,
                           '_links': {'self': _create_href(method_name='get_chromosome_assocs',
                                                           params={'chromosome': chromosome})}}
        chromosomes_list.append(chromosome_info)

    response = {'_embedded': {'chromosomes': chromosomes_list}}
    return simplejson.dumps(OrderedDict(response))


@app.route('/chromosomes/<string:chromosome>')
def get_chromosome_assocs(chromosome):
    args = request.args.to_dict()
    start, size, pval, pval_interval = _get_basic_arguments(args)
    study = _retrieve_endpoint_arguments(args, 'study_accession')
    bp = _retrieve_endpoint_arguments(args, 'bp')
    bp_interval = _get_interval(bp, IntInterval)

    searcher = search.Search(properties.output_path)

    try:
        datasets, index_marker = searcher.search_chromosome(chromosome=chromosome,
                                                            start=start, size=size, study=study,
                                                            pval_interval=pval_interval, bp_interval=bp_interval)

        data_dict = _get_array_to_display(datasets)
        params = {'chromosome': chromosome, 'p-value': pval, 'bp': bp, 'study_accession': study}
        response = _create_associations_response(method_name='get_chromosome_assocs', start=start, size=size, index_marker=index_marker,
                                                 data_dict=data_dict, params=params)

        return simplejson.dumps(OrderedDict(response), ignore_nan=True)

    except NotFoundError as error:
        raise RequestedNotFound(str(error))


@app.route('/variants')
@app.route('/variants/<string:variant>')
def get_variants(variant=None):
    if variant is None:
        raise InvalidUrl("Missing variant id")
    args = request.args.to_dict()
    start, size, pval, pval_interval = _get_basic_arguments(args)
    study = _retrieve_endpoint_arguments(args, "study_accession")
    chromosome = _retrieve_endpoint_arguments(args, "chromosome")
    searcher = search.Search(properties.output_path)
    if chromosome is None:
        raise ArgumentMissing(message="Required string parameter \'chromosome\' is missing")
    # try:
    datasets, index_marker = searcher.search_snp(snp=variant, chromosome=chromosome, start=start, size=size,
                                                 pval_interval=pval_interval, study=study)

    data_dict = _get_array_to_display(datasets, variant=variant)
    params = {'variant': variant, 'p-value': pval, 'study_accession': study}
    response = _create_associations_response(method_name='get_variants', start=start, size=size, index_marker=index_marker,
                                             data_dict=data_dict, params=params)

    return simplejson.dumps(OrderedDict(response), ignore_nan=True)

    # except ValueError:
    #     abort(404)


def main():
    _set_properties()
    app.run(host='0.0.0.0', port=8080)


if __name__ == '__main__':
    main()
