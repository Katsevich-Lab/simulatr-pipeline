/*
 * CHECK_CONTRACT — assert that the installed simulatr's contract
 * schema_version matches what the pipeline was built for.
 */

process CHECK_CONTRACT {
    tag 'contract'
    label 'CHECK_CONTRACT'

    output:
    path 'contract_check.txt'

    script:
    def expected = params.expected_contract_version ?: 1
    """
    check_contract.R ${expected} > contract_check.txt
    """
}
