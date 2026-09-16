# Research deployment in Germany

The interface supports German and English and German numeric formatting. This is independent research software, not a certified medical device, clinical decision-support system or Charité-endorsed platform. It has not undergone a formal GDPR, cybersecurity or medical-device conformity assessment.

## Implemented privacy measures

- New transplant analysis runs in browser memory. Its code contains no table-upload API calls, analytics SDK, external fonts or automatic persistence.
- Imports require coded identifiers and are validated before replacement. Full studies and reports are explicitly downloaded by the researcher.
- Synthetic data are clearly labeled. Importing external tables switches the study to imported/mixed; synthetic sources remain identified in individual records.
- A static build can be hosted on approved institutional infrastructure. Existing Python/Docker functionality remains available for the original omics workspace.

Local computation does not make genetic data anonymous. Source records, notes and downloaded reports may remain sensitive. Hosting access logs, authentication, IP addresses, browser extensions, device backups and downloaded files are outside this app's in-memory model. The managed demo's hosting location and contractual arrangements have not been verified as suitable for identifiable patient research data. Use the synthetic demo there; deploy internally and obtain institutional approval before real data.

## Institutional decisions before patient research

Discuss the specific study with the institution's data-protection officer and research governance team: controller/processor roles, GDPR Articles 6/9 legal basis, participant information and consent where applicable, ethics approvals, access roles, approved hosting location and contracts, transfer conditions, encryption, retention/deletion, incident procedures and whether a DPIA is required. Pseudonymization keys should be separated from analysis data. This prototype does not implement institutional RBAC, managed consent, tamper-evident audit logs or controlled retention.

[GDPR official text](https://eur-lex.europa.eu/eli/reg/2016/679/oj/eng) addresses special-category health/genetic data and research safeguards. [BDSG §27](https://www.gesetze-im-internet.de/bdsg_2018/__27.html) provides specific German research provisions with conditions and safeguards; it is not blanket permission for a project. Applicability must be assessed for the institution and study, including relevant state and hospital rules. No assertion of compliance is made by these engineering measures.

For translational/clinical use, establish intended purpose, regulatory classification, quality management, validation, traceability and professional oversight independently. A manual shortlist or functional-evidence label in this research prototype is not a clinical release decision.
